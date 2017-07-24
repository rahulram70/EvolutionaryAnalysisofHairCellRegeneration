# Check Orthology To DARP's
here we define what proteins share a 70% similarity to the zebrafish transcript
for each transcript.

## Methods
This experiment is broken down into 3 consecutive steps allowing for parts to be computed sequencially as time allows. this is valuable due to the length of time required to align our dataset.

| 1. homology_check.py |
|-------------------|
| this script generates alignment ready files conatining all proteins which are orthologous to a particular transcript | 
<br>

| 2. align.py |
|----------|
| this script takes the alignment ready files generated from the previous script and aligns those files, outputing to a seperate directory | 
<br>

| 3. comparison.py |
|---------------|
| this script takes the alignments from the previous script and compares all the proteins and returns a log of similarity between all proteins and the comparator species |
<br>

## Results

