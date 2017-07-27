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
    in_fl = res_dir + "/resources/data-cleaned/aligned_seq/ENSDART00000000042_0 .txt"

    out = pc.comp_cons_seq(in_fl)
    print(out)

    
if __name__ == '__main__':
    main()



