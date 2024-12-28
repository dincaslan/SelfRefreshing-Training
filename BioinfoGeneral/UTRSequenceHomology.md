If you ever used [NCBI-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), you might be familiar with "Sequence similarity (homology)" search.
I was also wondering whether there is a way to compare 5 prime UTR sequence similarity of all transcripts of ortholog* genes of mouse and human (you can expand to other species as well.)
If you are unsure about orthology, then plese feel free to use [Emsemb BiomaRt](https://www.ensembl.org/info/data/biomart/index.html) to check it out (it has both an online tool and bioconductor package to run in R).

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
blastn -query /Users/dincaslan/Desktop/HS_FLT1.txt -subject /Users/dincaslan/Desktop/MS_FLT1.txt -out /Users/dincaslan/Desktop/test_blastn.txt 
```

```r
```


```r
```


```r
```


```r
```


```r
```
