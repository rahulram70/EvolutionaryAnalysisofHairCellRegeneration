#!/usr/bin/env python3

from __future__ import absolute_import
import os
import procomp as pc
import time

@pc.stats
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    inputFile = script_dir + "/procomp-out.txt"
    outputFile = script_dir + "/procomp-analyzed.txt"


    #pc.MainPC_analysis(inputFile, outputFile)
    

if __name__ == '__main__':
    main()

    
