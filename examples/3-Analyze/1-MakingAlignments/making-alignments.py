#!/usr/bin/env python3

from __future__ import absolute_import
import os
import procomp as pc
import time

@pc.stats
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    aligned_outF = script_dir + "/alignedSeqOutput/"
    unalignedF = script_dir + "/unalignedSeqs/"

    mspath = os.path.abspath(os.path.join(script_dir, os.pardir))
    mspath = os.path.abspath(os.path.join(mspath, os.pardir))
    mspath += "/resources/muscle3.8.31_i86darwin64/muscle3.8.31_i86darwin64.31_i86darwin64"
    
    pc.bioMuscleAlign(unalignedF, mspath, aligned_outF)

if __name__ == '__main__':
    main()

    
