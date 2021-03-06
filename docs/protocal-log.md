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
            - rows now: 11339, rows starting: 12711
            - 1372 rows removed

### July 19, 2017
- Project Meeting discussed solutions for the 1-to-many issue regarding orthology due to the fact that combinations generated for roughly 10% of the transcript result in un-computably large numbers of combinations (reaching >1x10^60).
    - Solution 1
        - before computing combinations, align regenerators and non-regenerators seperately to identify a conserved sequence for each group by finding proteins with >70% similarity in the filter alignment stage.
    - Solution 2
        - if the alignment filter stage does not significantly reduce the number of combinations for a specific transcript further refactoring methods must be considered which will be discussed at a later point.
        
### August 1, 2017
- Project meeting to discuss how to handle combinations and analysis. Robby, Alli, Phil, Rahul, and Matt Lambert in attendance.
- Need to measure how many genes we are dealing with. Take transcript list and convert to gene ids
- Ways to widdle down combinations
    - Preferred way is to go with the longest sequence for each species and transcript header making it a one to one
    - 40% sequence homology is too stringent for vertebrates 
    - Determine whats a truncated protein to elminate some noise?
- Recreating ancestral sequences for analysis 
    - Create ancestral sequences for fish, (birds + reptiles), and mammals using PAML. 
        - For mammals consider excluding monotremes
- Use Amazon webservices for PAML or other process intensive computations
- Consider analyzing gene duplicates that are conserved in the regenerator but different or lost in the non-regenerators.

### August 2, 2017
- Project meeting to discuss ortholog approach. Robby, Matt, and Phil in attendance.
- New strategy to obtain sequences and determine orthologs.
    - Use NCBI's database (allows for reciprocal blast determination of orthologs and more diverse species).
    - Obtain zebrafish sequences associated with every gene from RNA seq papers.
    - Keep only longest zebrafish protein sequence for each gene id.
    - Blast those sequences (blastp) and obtain best hit for species list that Phil will select.
        - Choose well annotated species with broad distrubution across species phylogeny (ideally scaffolds=chromosomes).
        - Should there be a threshold for the blast (Either the first or the reciprocal)?
    - Then take those sequences and blast them back to zebrafish.
        - If the blast back to zebrafish does not yield the orginal query sequence then no ture ortholog exists for that pair.
        - This result will yield at most 1 sequence per species per gene id.
- Discussed analysis methods.
    - Analyze avergae number of acceptable amino acid substitutions per site (regens vs. non).
    - Stringent positive selection vs. relaxation of seleciton.
    - PAML analysis.
        - Control file is analogous to config file.
        - Need to create a gene tree for each gene.
        - Not great to analyze alot of genes because the nodes where inferred ancestral sequence are need to be manaually selected.
