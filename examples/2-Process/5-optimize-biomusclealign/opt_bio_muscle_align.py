#!/usr/bin/env python3

import os
import sys

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    # go up 2 directories to project dir
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))

    # add project to path and import package
    sys.path.append(res_dir)
    import procomp as pc

    # Setup paths to needed directories
    #
    
    out_dir = script_dir + "/post-align/"
    in_dir = script_dir + "/pre-align/"
    musl_path = res_dir + "/resources/muscle/muscle3.8.31_i86darwin64.31_i86darwin64"
    
    # ---------------------------------
    # original muscle align function
    #
    pc.bioMuscleAlign(in_dir, musl_path, outputF=out_dir)

    # ------------------------------
    # New multithreading based alignments
    #
    


    
if __name__ == '__main__':
    main()