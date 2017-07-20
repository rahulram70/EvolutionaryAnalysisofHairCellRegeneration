# Removing Transcript IDs that have either have 0 non-regenerating or regenerating species or that are 1 to 1 comparisons
Here we remove all instances of transcript IDs that have either have 0 non-regenerating or regenerating species or that are 1 to 1 comparisons. However, we do not edit the ex_out_comb text file but add the edits into a new file called revised_ex_out_comb text file.

## Methods
In order to remove the Transcript IDs, we create a path to the ex_out_comb text file and the revised_ex_out_comb text file. Then we loop through the ex_out_comb file and identify the transcript headers. Then we analyze the headers to determine that it doesn't contain 0 non-regenerating species or
0 regenerating species. We also analyze it to determine if it is not a 1 to 1 comparison. Once those criteria are met the list is written into a separate text file called revised_ex_out_comb with its 
subsequent proteins.

## Results
The cleaned transcript IDs can be viewed in the revised_ex_comb text file in the resources folder.