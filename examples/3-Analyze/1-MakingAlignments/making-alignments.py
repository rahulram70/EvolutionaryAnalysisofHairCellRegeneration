#!/usr/bin/env python3

from __future__ import absolute_import
import os
import procomp as pc

@pc.stats
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    # go up 2 directories
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))


    out_dir = res_dir + "/resources/data-cleaned/aligned_seq/"
    in_dir = res_dir + "/resources/data-cleaned/unaligned_seq/"

    mspath = res_dir + "/resources/muscle3.8.31_i86darwin64/muscle3.8.31_i86darwin64.31_i86darwin64"
    

    pc.bioMuscleAlign(in_dir, mspath, out_dir)

if __name__ == '__main__':
    main()

    
