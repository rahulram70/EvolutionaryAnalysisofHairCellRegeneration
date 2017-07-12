# !/usr/bin/env python3

from __future__ import absolute_import
import os
import procomp as pc
import process as ps
import time

@pc.stats
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    inputFile = os.path.abspath(os.path.join(script_dir, os.pardir))
    inputFile = os.path.abspath(os.path.join(inputFile, os.pardir))
    outputFile = inputFile + "/resources/Jiang_2014_Raw-out.csv"
    inputFile += "/resources/Jiang 2014 Raw.csv"
    
    ps.row_rm_by_col_cond (inputFile, col='flag', fl_out=outoutputFile, cond="0")
    
    

if __name__ == '__main__':
    main()

