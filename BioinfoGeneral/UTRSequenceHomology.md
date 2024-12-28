If you ever used [NCBI-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), you might be familiar with "Sequence similarity (homology)" search.

### Using Command Line
The first way is using command lines: wget, eutils and blast+. 
Since most tools has whole FASTA and the coding sequence only FASTA (CDS), I simply substract the CDS from Full FASTA to have 5 prime non-coding sequences.

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

Another good news is that you can actullay use ncbi-blast tool in the command line (terrific!).
I used [homebrew](https://formulae.brew.sh/formula/blast) to downlaod the Blast+.

```bash
# Let's say you want to compare human and mouse ACTB gene (you can adjust according to your files)
# The file of the format should be FASTA for blast+ to work
# ">NM_002019.4 Homo sapiens fms related receptor tyrosine kinase 1 (FLT1), transcript variant 1, mRNA
# ATCGAGGTCCGCGGGAGGCTCGGAGCGCGCCAGGCGGACACTCCTCTCGGCTCCTCCCCGGCAGCGGCGG....

blastn -query /Your/Path/To/The/File/HS_ACTB.txt -subject /Your/Path/To/The/File/MS_ACTB.txt -out /Your/Path/To/The/File/test_blastn.txt

# The output format is not a single score. You need to think about how to extract the relevant information in a single file (I skip this since we have an alternative below).
```

Here is a useful [BLAST+ tutorial](https://conmeehan.github.io/blast+tutorial.html).

### Using R
I was also wondering whether there is a way to compare 5 prime UTR sequence similarity of all transcripts of ortholog* genes of mouse and human (you can expand to other species as well.)
The alternative is using Ensembl-BiomaRt for the retrieval of the relevant non-coding/coding FASTA and compare them using stringdist package of R for the similarity.
> If you are unsure about orthology, then plese feel free to use [Emsemb BiomaRt](https://www.ensembl.org/info/data/biomart/index.html) to check it out (it has both an online tool and bioconductor package to run in R).
```r
# Downloading the relevant packages
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
library(biomaRt)

### Thanks ChatGPT for reminding me BiomaRt!, and here is nice biconductor tutorial if you are curious: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html ###
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
mtranscripts$mgi_symbol<-toupper(mtranscripts$mgi_symbol) # make mgi symbol similar to hgnc for downstream processes
names(mtranscripts)[names(mtranscripts) == "mgi_symbol"] <- "hgnc_symbol"

# Initialize a data frame to keep the results
results <- data.frame(
  hGoI = character(),
  hToI = character(),
  mGoI = character(),
  mToI = character(),
  hgnc_symbol = character(),
  Distance = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each element in gene of interest of two selected species to calculate distance, and save it
for (genes in unique(htranscripts$hgnc_symbol)) {
  # Filter X and Y for the current group
  hGoI <- htranscripts[htranscripts$hgnc_symbol == genes,]
  mGoI <- mtranscripts[mtranscripts$hgnc_symbol == genes,]
  
  # Proceed only if there are matches in both X and Y for this group
  if (nrow(hGoI) > 0 && nrow(mGoI) > 0) {
    # Compare each sequence in X with each sequence in Y
    for (i in 1:nrow(hGoI)) {
      for (j in 1:nrow(mGoI)) {
        # Extract sequences and types
        h_type <- hGoI$ensembl_transcript_id[i]
        h_seq <- hGoI[,"5utr"][i]
        m_type <- mGoI$ensembl_transcript_id[j]
        m_seq <- mGoI[,"5utr"][j]
        
        dist <- stringdist::stringdist(h_seq, m_seq, method = "lv") # "lv" for Levenshtein distance:
# "counting the weighted number of insertions, deletions and substitutions necessary to turn one string into another"
        
        # Append the results data with each comparsion
        results <- rbind(results, data.frame(
          hGoI = h_type,
          hToI = h_seq,
          mGoI = m_type,
          mToI = m_seq,
          hgnc_symbol = genes,
          Distance = dist
        )) ## It would be more meaningful if there is a treshold for the distance measurements based on the mismatches.
      }
    }
  }
}

head(results) # only first 6 row of final results 
```
Here I attached the [Quarto doc]() of how it looks like. 
