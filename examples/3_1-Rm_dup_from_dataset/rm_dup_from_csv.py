#!/usr/bin/env python3

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

    outputFile = inputFile + "/resources/Jiang_and_Steiner Merged Transcript list-out.csv"
    inputFile += "/resources/Jiang and Steiner Merged Transcript list.csv"
    
    ps.row_rm_by_dup (inputFile, 'Transcript stable ID', outputFile)

if __name__ == '__main__':
    main()

    
