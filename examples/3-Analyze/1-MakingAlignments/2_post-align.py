#!/usr/bin/env python3

from __future__ import absolute_import
import os
import sys
import procomp as pc

@pc.stats
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    # go up 2 directories
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))

    L_in = sys.argv
    for i in range( len(L_in) ):
        if L_in[i] == "-in":
            try:
                in_dir = L_in[i+1]
            except:
                print("ERROR: no input for -in")
        elif L_in[i] == "-out":
            try:
                out_dir = L_in[i+1]
            except:
                print("ERROR: no input for -out" )
        elif L_in[i] == "-mus":
            try:
                mspath = L_in[i+1]
            except:
                print( "ERROR: no input for -mus" )

    """
    # the folder where the alignments will be written too
    out_dir = res_dir + "/resources/data-cleaned/aligned_seq/"
    

    # the folder where the pre-alignments are located
    in_dir = res_dir + "/resources/data-cleaned/unaligned_seq/"

    # the path for the muscle binary file
    mspath = res_dir + "/resources/muscle3.8.31_i86darwin64/muscle3.8.31_i86darwin64.31_i86darwin64"
    """

    # compute the alignment and output results to output folder
    pc.bioMuscleAlign(in_dir, mspath, out_dir)

if __name__ == '__main__':
    main()

    
