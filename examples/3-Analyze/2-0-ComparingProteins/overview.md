## Using Procomp to analyze the difference in aligned sequences between regenerators and non-regenerators

Here we demonstrate a central analysis to this project. aligned amino sequences are analyzed for the percentage of similarity between species in the same category (defined by the species-groups.txt) and difference between categories. The idea is that proteins containing a large difference between the two groups in our study may indicate a candidate which preserved function for regernation while the emerging mammalian branch evolved away from that function. 

## Methods

|compare_pro_in_folder.py|
|------------------------|
|our method for comparing aligned proteins follows a 3 step procedure. |    
|1. foreach alignment, iterate over every site and add all sites from group 1 to a list and all sites from group 2 to a list |
|2. compare lists for each group, if the groups fit a number of criteron add 1 to the hit count |
|3. output the total number of hits and information about the protein |


## Results

Output currently follows the format:

|  site = 419  &nbsp;&nbsp;&nbsp;&nbsp; R: EEEEE &nbsp;&nbsp;&nbsp;&nbsp; NR: AAADD<br>
|  site = 471  &nbsp;&nbsp;&nbsp;&nbsp; R: AAAAA &nbsp;&nbsp;&nbsp;&nbsp; NR: CYCCC<br>
|  site = 477  &nbsp;&nbsp;&nbsp;&nbsp; R: HHHHH &nbsp;&nbsp;&nbsp;&nbsp; NR: EEEEE<br>
L__ GENE: ENSDARG00000a44125_10.txt  &nbsp;&nbsp;&nbsp;&nbsp;---total Length:539 &nbsp;&nbsp;&nbsp;&nbsp;---#hits =3 &nbsp;&nbsp;&nbsp;&nbsp;R: 5 &nbsp;&nbsp;&nbsp;&nbsp;NR: 5<br>