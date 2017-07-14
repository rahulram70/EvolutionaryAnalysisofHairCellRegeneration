## Comparing the Ensembl Database (ED) to the Current-transcript IDs (CTI)

 This program compares the protein IDs fron the Current Transcripts to the Ensembl Databse. This program checks
 to see if an ID from the CTI is located in the Ensembl Database. It also checks if the ID in the ED has a correspoinding protein sequence that ends in an *. If there's no ID or a protein sequence error a list is returned to the user containing the erroneous protein IDs.

 ## Methods

The program breaks down the comparison into 3 parts. 

1. Retrieve the CTI, places them into a list, and return it to the user
2. Retrives the ED, separate the IDs and Sequences and place them into lists, and return it to the user
3. Loop through the CTI ID's and determine if the ID is found in the ED. Also, check to see if each ID in the ED contains a correct protein sequence. For those ID's that do not match both requirements, they are     placed in an error list and returned to the user.

** For a more indepth understanding, please look at comparison.py and read the comments

## Results

    Program was successful, below are the more detalied resutls:

    List of Protein ID's with incorrect sequences: 
    
    ['ENSPFOP00000009023', 'ENSPFOP00000009023', 'ENSPFOP00000009023', 'ENSPFOP00000009023', 'ENSGMOP00000008763', 'ENSGMOP00000008763', 'ENSGMOP00000008763', 'ENSMPUP00000009176', 'ENSMMUP00000055636', 'ENSMMUP00000055636', 'ENSMMUP00000055636', 'ENSMMUP00000055636', 'ENSRNOP00000070247', 'ENSRNOP00000070247', 'ENSGACP00000001929', 'ENSGACP00000001929', 'ENSGACP00000001929']

