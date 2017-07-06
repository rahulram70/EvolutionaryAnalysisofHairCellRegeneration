## Transcript id checks

Here we analyze the transcript ids gathered from ensembl for each species to verify that no transcript are missing from a particular speices that would inaccurately influence the output of calculating combinations and subsequently performing our main amino acid analysis on the dataset.

## Methods

|tr_id_check_mstr.py|
|---|
|this script is used to print out results from analyzing the master id list <i>(/resources/data-cleaned/Jiang_and_Steiner Merged Transcript list-out.csv)</i>. outputs from this unit test are stored in seperate text files.|


## Results

|file_and_id_check.txt|
|---|
|the test verified that for each species we had equal number of files and transcript ids. the master id list contains 11,339 transcript ids and all species contain 10,170 ids. there are 1,169 ids which do not show up in the local lists for each species due to the fact that ensembl does not have any information for these transcripts. this list is maintained in the accompanying file <i>(ids_missing_in_all.txt)</i>
<br>

|ids_missing_in_all.txt|
|---|
|the ids in this list are ids present in the master transcript id list but not present in the local list of ids for all the species. further investigation found this list of ids to be empty entries on ensembl.|
