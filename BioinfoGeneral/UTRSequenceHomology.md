If you ever used [NCBI-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), you might be familiar with "Sequence similarity (homology)" search.

### Using Command Line
```bash
# More info about eutils and efectch, https://www.ncbi.nlm.nih.gov/books/NBK179288/, https://github.com/NCBI-Hackathons/EDirectCookbook, https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/, https://www.ncbi.nlm.nih.gov/books/NBK25499/,
# and ncbi datasets, https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/, https://github.com/ncbi/datasets, https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/gene-package/, https://www.nature.com/articles/s41597-024-03571-y  as an alternative

# "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=[NM TRANSCRIPT ID]/&rettype=[FASTA FORMAT]"
URL_base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=%s"

CompleteFASTA=$(wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_010240.2&rettype=fasta")
CDSonly=$(wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_010240.2&rettype=fasta_cds_na")

echo "$CompleteFASTA" | grep ">" | while read -r header; do
	accession=$(echo "$header" | sed -E 's/^>[^|]*\|([^ ]+).*/\1/')
	location=$(echo "$header" | sed -n 's/.*location=\([0-9]*\)\.\..*/\1/p')
  # Identify the CDS region from whole sequence and leave the non-coding (nc)
	seq_in_whole=$(echo "$CompleteFASTA" | grep -A 1 "$accession" | tail -n 1)
	start_pos=$(echo "$location" | cut -d. -f1)
	UTR5p_seq="${seq_in_whole:0:$((start_pos-1))}" # leaves the 5UTR nc sequnece

	echo ">${accession}_UTR5p" >> $output_file # not necessary but keep in similar accession, gene name etc. format
	echo "$UTR5p_seq" >> $output_file
	echo "Completed: $accession to $output_file"
done

```
```bash

# You can try a for llop for multiple inputs 

# out_file="UTR5p_out.txt"
# TranscriptIDs=("NM_010240.1", "NM_010240.2") # examples
# FASTAform=("fasta", "fasta_cds_na") # to substract cds from the complete sequence to have 5 and 3UTR depending on interests, mine is 5pUTR

# for id in "${TranscriptIDs[@]}"; do
#   for fform in "${FASTAform[@]}"; do
#     URL=$(printf "$URL_base" "$id" "$fform")
#     wget "$URL" -O ->> $out_file
#   done
# done
```
Another good news is that you can actullay use ncbi-blast tool in the command line (terrific!).
I used [homebrew](https://formulae.brew.sh/formula/blast) to downlaod the Blast+
[BLAST+ tutorial](https://conmeehan.github.io/blast+tutorial.html)
```bash
# Let's say you want to compare human and mouse ACTB gene
blastn -query /Your/Path/To/The/File/HS_ACTB.txt -subject /Your/Path/To/The/File/MS_ACTB.txt -out /Your/Path/To/The/File/test_blastn.txt 
```
### Using R
The alternative is using Ensembl-BiomaRt for the retrieval of the relevant non-coding/coding FASTA and compare them using stringdist package of R for the similarity.
I was also wondering whether there is a way to compare 5 prime UTR sequence similarity of all transcripts of ortholog* genes of mouse and human (you can expand to other species as well.)
If you are unsure about orthology, then plese feel free to use [Emsemb BiomaRt](https://www.ensembl.org/info/data/biomart/index.html) to check it out (it has both an online tool and bioconductor package to run in R).

```r
# Downloading the relevant packages
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
library(biomaRt)

### Thanks ChatGPT for reminding me BiomaRt! ###
# Connect to EnsemblDB and select dataset genes (human and mouse)
hensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
```


```r
# List genes of interests
hGoF <- c('FLT1', 'ACTB') # human
mGoF <- c('Flt1', 'Actb') # mouse

# Retrieve transcript accessions for GoI
htranscripts <- getBM(
  attributes = c("hgnc_symbol", "ensembl_transcript_id", "5utr"),
  filters = "hgnc_symbol",
  values = hGoF,
  mart = hensembl
)

htranscripts # see the list

mtranscripts <- getBM(
  attributes = c("mgi_symbol", "ensembl_transcript_id", "5utr"),
  filters = "mgi_symbol",
  values = mGoF,
  mart = mensembl
)

mtranscripts # see the list

# Remove the colums with fasta sequences unavailable
htranscripts <-htranscripts[ htranscripts["5utr"]!="Sequence unavailable",]
mtranscripts <-mtranscripts[ mtranscripts["5utr"]!="Sequence unavailable",]

# Save it as a CSV file if you wish
# write.csv(htranscript_data, "htranscript_accessions.csv", row.names = FALSE)
# write.csv(mtranscript_data, "mtranscript_accessions.csv", row.names = FALSE)
```


```r
# Approximate string matching, the similarity between two strings (the paper: https://cran.r-project.org/web/packages/stringdist/vignettes/RJournal_6_111-122-2014.pdf)
# Initialize a data frame to keep the results
results <- data.frame(Human = character(), Mouse = character(), Distance = numeric(), stringsAsFactors = FALSE)

# Loop over each element in gene of interest of two selected species to calculate distance, and save it
for (x in htranscripts$ensembl_transcript_id) {
  for (y in mtranscripts$ensembl_transcript_id) {
    dist <- stringdist::stringdist(x, y, method = "lv") # "lv" for Levenshtein distance: "counting the weighted number of insertions, deletions and substitutions necessary to turn one string into another"
    results <- rbind(results, data.frame(Human = x, Mouse = y, Distance = dist))
  }
}

# It would be more meaningful if there is a treshold for the distance measurements based on the mismatches
results # final results
```
