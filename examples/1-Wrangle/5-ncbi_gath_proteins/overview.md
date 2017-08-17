# Gathering Zebrafish Protein Sequences from NCBI
This example focusses on utilizing the Entrez api provided by NCBI via 
the wrapper in biopython, to download all the protein sequences for the
zebrafish genes.

## Methods

|get_gene_names.pl |
|-|
|this script uses the ensembl perl api to aquire the external gene names from the zebrafish ensembl gene id |
<br>

|gath_proteins.py |
|-|
|this script uses the biopython wrapper for the ncbi Entrez api to query for the gene names collected from get_gene_names.pl and then download their protein translations as well as info related information about the proteins |
<br>

|reciprocal_blast.py |
|-|
|after zebrafish proteins have been downloaded we use this script to run a reciprocal blast on every zebrafish protein and download all proteins which have the original zebrafish protein as a top hit from the blast |
<br>

## Results
