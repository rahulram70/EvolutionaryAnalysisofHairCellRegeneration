#!/usr/bin/env python3

import os
import sys

print(__file__)

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
    unalign_path = res_dir + "/examples/2-Process/5-optimize-biomusclealign/pre-align/"
    align_path = res_dir + "/examples/2-Process/5-optimize-biomusclealign/post-align/"
    muscle_path = res_dir + "/resources/muscle/muscle3.8.31_i86darwin64.31_i86darwin64"
    

    pc.bioMuscleAlign(align_path, muscle_path)


if __name__ == '__main__':
    main()