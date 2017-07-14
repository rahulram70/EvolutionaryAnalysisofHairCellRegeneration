
# Removing Duplicate Proteins from Downloaded Data
Here we remove all instances of duplicate proteins and their sequences from the data downloaded from ensembl in the protein_sequences directory. This is neccessary because proteins with same id's also have same sequences which will result in longest search queries when the combinations files are created.

Duplicate Proteins exist because some proteins have are orthologous to multiple transcripts. Because of this iterating through all the protein ids disregarding its ortholog transcript id resulted downloading some proteins more than once unfortunately.

## Methods
in order to remove the duplicate proteins we first open each species file in /resources/data-raw/protein-sequences/ and linearly iterate through it adding proteins to a list if they do not exist in that list, then write that list to a file when all proteins have been checked. this consequently excludes any duplicate proteins that may have existed in the originally downloaded data.

## Results
The cleaned protein sequences can be viewed at<br> <i>/resources/data-cleaned/protein-seq/</i>