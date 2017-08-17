#!/usr/bin/env python3

import os
import sys
import time

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
    
    filePath = res_dir + "/resources/data-raw NCBI/NCBI_Species_list.xlsx"

    print(pc.gen_dataframe(filePath))    
if __name__ == '__main__':
    main()