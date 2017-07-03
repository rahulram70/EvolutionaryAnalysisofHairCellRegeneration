## Clean dataset of genes by removing insignificant and duplicate entries

Here we take a raw list of genes from the paper (Jiang 2014) and eliminate all entries containing zeros in the flag (standard deviation thresholds) column. the processed dataset can then be used to aquire transcript id's which will ultimately be used to generate combinations of ortholog proteins.

<br>

## Methods

|clean_dataset.py |
|-----------------|
|here we input a .csv file and convert to a pandas.dataframe object then iterate through all rows, eliminating entries containing the conditional (cond). the conditional for this portion was the flag column which contained standard deviation thresholds. |
<br>

|rm_dup_from_csv.py |
|-------------------|
|this script takes a .csv file and returns a .csv file which has duplicate values in a particular column removed. here we remove duplicates from the 'transcript stable id' column. |
<br>


## Results

processing is successful and results in a .csv file in the resources directory.