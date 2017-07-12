## Creating combinations of ortholog proteins

Here we generate combinations of orthologous proteins that correlate to Zebrafish transcripts and output a file of combinations in a semi-readable format.

## Methods

|creating-combinations.py|
|---|
|this script utilizes the 3 core functions in procomp for generating combinations of protein ids.|
|<h4>comb_TrPr()</h4> first we combine all transcript or genes and the proteins they relate to, into a single list.|
|<h4>comb_rm_dups()</h4> the list is then analyzed, removing duplicate protein instances as well as generating a makeshift table containing all ortholog proteins for a given transcript or gene id. this is done without calculating combinations.| 
|<h4>comb_gen_combs()</h4> here, the ortholog combinations count are first calculated. if below a set amount, the combinations are written to an output file, otherwise the species with the most orthologous proteins for the given gene or transcript is eliminated from the list, then the combinations are recalculated. this refactoring continues until the the combination count is below a given threshold.|

#### REFACTORING TRANSCRIPTS
There occured a problem with some transcript which generated too many combinations due to the fact that large numbers of orthologous proteins. In order to address this issue, and bring the combination count down to a workable count we use a refactor method which takes transcripts with combinations over 500 and eliminates species

## Results

Output is written to a txt file in the format:<br>
|'transcript id' | 'species 1 protein id' | 'species 2 protein id'| ect. |
| - | - | - | - |
|ENSDART00000000221 |ENSACAP00000015676 |ENSAMEP00000005565 | ... |
|'                ' |ENSACAP00000015769 |ENSAMEP00000005592 | ... |
|'                ' |ENSACAP00000017586 |ENSAMEP00000005646 | ... |