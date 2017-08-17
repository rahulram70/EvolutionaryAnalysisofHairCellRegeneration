# Local Blast Setup and Execution
 
## Contents 

- [Dependancies](#Dependancies)
- [Setting Up](#setting-up)

<hr>

## Dependancies

In order to perform a local blast 2 key elements are needed.
- [NCBI Blast+ Suite](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Database of Sequences](https://www.ncbi.nlm.nih.gov/protein)

In our case we created our own database of protein sequences by taking all protein sequences for our list of species.
<hr>

## Usage

First the Blast+ Suite must be installed as well as at least two sequence files. after these resources have been collected a database must be setup for each species.

### 1. SETUP DATABASE
``` 
makeblastdb -dbtype prot -in gallus_gallus.fasta -out Gallus__db -parse_seqids 
```
|variable | description |
|-|-|
|gallus_gallus.fasta |the file with all the fasta entries for the speices |
|Gallus__db |what to name the database |
<br>

### 2. RUN BLASTP
```
blastp -query ZB.fasta -db Gallus__db -out ZF_vs_gallus.bls -num_alignments 1 -num_descriptions 1
```
|variable | description |
|-|-|
|ZB.fasta | the file with all input species fasta entries |
|Gallus__db |the species database to query against |
|ZF_vs_gallus.bls |the output file path for the blast results |
|-num_alignments |number of hits to return from blast |
<br>

    