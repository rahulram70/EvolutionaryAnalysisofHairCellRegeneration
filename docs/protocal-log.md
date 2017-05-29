# Protocol Log
This file keeps track of the protocals used and logic
behind the project progression.

### May 25, 2017
- Pulled differentially regulated transcripts from two PNAS papers (Jiang et al 2014 and Steiner et al 2014). 
    - From Jiang paper supplemental excel file 2. All “flagged” genes (-5,-3,-1,1,3,5) were used. Robby excluded 0’s from the list with code.
        - In resources folder on procomp “Jiang_2014_Raw-cleaned
        - From there input 500 gene id’s at a time into ensembl biomart to get associated transcript id’s
    - From Steiner paper Dataset S4 was used in entirety
    - Merged two transcript id lists
        - Robby removed duplicates

