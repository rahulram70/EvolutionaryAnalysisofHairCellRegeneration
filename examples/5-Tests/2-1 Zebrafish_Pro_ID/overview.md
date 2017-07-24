## Checking for Missing Protein IDs
The purpose of this script was to loop through the entire mst_D dictionary and return the transcript IDs that do not have a zebrafish protein ID or had multiple zebrafish protein IDs.

## Methods
We performed this action by opening the Species Group file, converting it into a list, and removing the duplicates from it. Then we ran the list in the duplicate remover function in procomp to create a table and remove duplicate proteins. Next, the generated table was converted into a dictionary with the keys being defined as the zebrafish transcript IDs, and the values being defined as its subsequent protein IDs. Finally, we looped through the dictionary and returned transcript IDs that did not have a zebrafish protein ID or contained multiple zebrafish protein IDs.

## Results
After running the program, we have determined that 1086 zebrafish transcript IDs do not contain proteins. Also, we have determined that two transcript IDs contain multiple proteins. They are ENSDART00000115071 and ENSDART00000125016. These IDs are duplicates with the same transcript and protein IDs.