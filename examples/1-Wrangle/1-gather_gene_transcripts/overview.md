## Gathering transcripts and protein sequences from gene ids in Ensembl

This experiment aims to retrieve protein sequences given the transcript or gene id. this task was approached using, 1) the ensembl perl api and, 2) the python based genomic tools library. These script will allow for faster data retrieval from the ensembl database resulting in much time saved.

<br>

## Methods

| ens_tr_to_pr.pl |
|-----------------|
|this script cycles through the directories containing raw transcript-protein associations and uses the ensembl core database to query the database for the given transcript and then retrieve the translation object thereby aquiring the desired protein id and the unalinged protein sequence. |
<br>

| gather_gene_tr.py |
|-----------------|
|this script uses pycogent to query ensembl for a gene object and then print out all exons associated with each transcript of the gene. |
|NOTE: this was for practice implementations only and may not be included in the final publication of this project |
<br>

| test.pl |
|-----------------|
|this script uses the ensembl compara database to query for a protein (family adapter object) and then also retrieve all ortholog proteins relating to that protein. |
|NOTE: this was for practice implementations only and may not be included in the final publication of this project |
<br>
<br>

## Results

this project is still underway... a throughput implementation of either pycogent or the ensembl perl api has not been developped yet due to the ambiguity of the ensembl api.