
__name__ = "Process Module"

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

def row_rm_by_col_cond (fl_in, col, fl_out="", cond="0"):
    """
    Overview:
        Removes rows from a dataset by the condition specified
        in given column.        
    Inputs:
        in_fl   = file path of csv to be processed
        out_fl  = file path of csv to be written to
        col     = title of the column to be checked
        cond    = condition required for the row to be 
                  removed. default set to "0".
    Outputs:
        modified csv will output to a new file called 'X-out.csv'
        unless no output file is specified, then output is written
        to the input file.
    """

    if (fl_out == ""):
        fl_out = fl_in

    df  = pd.read_csv(fl_in)
    col = df.columns.get_loc(col)    
    L   = []
    
    # 1. make a list of which rows to drop
    for index, row in df.iterrows():
        val = df.iloc[index, col]
        val = val.split(",")
        if (cond in val):
            L.append(index)

    # 2. removing elements from L 
    df = df.drop(L)

    # 3. convert dataframe to csv output file
    df.to_csv(fl_out)

def row_rm_by_dup (fl_in, col, fl_out=""):
    """
    Overview:
        Removes rows from data frame that contain duplicate values
        in a particular column.
    Inputs:
        f_in    = csv file path
        col     = title of the column to be checked
        f_out   = csv file path to write results to.
    Output:
        modified csv will output to a new file called 'X-out.csv'
        unless no output file is specified, then output is written
        to the input file.
    """
     
    if (fl_out == ""):
        fl_out = f_in
    df  = pd.read_csv(fl_in) 
    initLen = df.shape[0]

    # 1. sort dataframe by selected column
    df.sort_values(col, inplace=True)

    # 2. remove elements from L
    df = df.drop_duplicates()
 
    # 3. reindex list since values have been removed
    df = df.reset_index()
    df = df.drop("index", axis=1)

    # 4. convert dataframe to csv output file
    print("rows now: {}, rows starting: {}".format(df.shape[0], initLen))
    print("{} rows removed".format( (int(initLen) - int(df.shape[0])) ))
    df.to_csv(fl_out)

def merge_folder_to_file(fd_in, f_out, save=1):
    """
    Overview
        Compile a folder of text files into a single text file.
    Inputs:
        fd_in   = folder of text files
        f_out   = file to compile files to
        save    = option to save old text files
    Outputs:
        f_out   = ideally put the output file outside the folder
                  being compiled from.
    """

    L = []
    for file in os.listdir(fd_in):
        if file.endswith(".txt"):
            path = fd_in + file
            val = open(path, "r").read()
            L.append(val)

    with open(f_out, "w+") as out:
        for i in L:
            out.write(i + '\n')

