#!/usr/bin/env python3

from __future__ import absolute_import
import os
import procomp as pc
import time

@pc.stats
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    alignedFolder = script_dir + "/demo-alignments/"
    outputTxt = script_dir + "/procomp-out.txt"
    ids = script_dir + "/sp_id-Jan6-2017.txt"
    
    pc.MainProteinCompare(alignedFolder, ids, "R", "NR", outputTxt)
    

if __name__ == '__main__':
    main()

    
