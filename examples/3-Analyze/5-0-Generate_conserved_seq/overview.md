# Defining Orthologs by aligning groups for each transcript
The goal of this example is to reduce the number of combinations generated from
each transcript by first aligning all proteins for each group (R and NR) then
determining which proteins meet a 70% similarity match to represent a conserved
sequence. proteins not matching the 70% similarity will be excluded from analysis
and reduce the number of combinations of orthologs.

## Methods
This task can be broken down into multiple parts

1. make 2 lists, R and NR
2. align both groups seperately
3. check each protein in each alignment for a 70% similarity to the other species 

## Results

