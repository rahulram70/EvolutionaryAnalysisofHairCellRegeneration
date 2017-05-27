
"""
# Author:      Robby Boney, robby.boney@wsu.edu
# License:     ...
# Last Update: May 26, 2017
# Start Date:  May 26, 2017
# Version:     alpha 1.1

This script contains functions for cleaning up data,
and processing inputs & outputs to be used with 
procomp functions.
"""

import pandas as pd
import os
from os import listdir
from os.path import isfile, join

def row_rm_by_col_cond (f_in, col, cond, f_out=""):
    """
    Removes rows from a dataset by the condition specified
    in given column.
    
    Inputs:
        f    = csv file path
        col  = title of the column to be checked
        cond = condition required for the row to be 
               removed.
    Outputs:
        modified csv will output to a new file called 'X-out.csv'
    """

    if (f_out == ""):
        f_out = f_in

    df = pd.read_csv(f_in)
    #print(df)
    col = df.columns.get_loc(col)
    
    L = []
    
    print("1. starting")
    for index, row in df.iterrows():
        val = df.iloc[index, col]
        val = val.split(",")
        if ('0' in val):
            L.append(index)

    print("2. removing {} elements".format(len(L)))
    df = df.drop(L)
    
    #i = 0
    #while i < len(L):
        #print(df.index[L[i]])
    #    )        
    #    i += 1

    print("3. converting {}".format(df.shape))
    
    df.to_csv(f_out)
    return 0

    
