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
    
    filePath = res_dir + "/resources/data-raw NCBI/BlastResults-ZF-to-X/ZF_vs_gallus.bls"
    output = res_dir + "/resources/data-raw NCBI/pro-seq/gallus gallus.fasta"
    inputL = []
    # ---------------------------------
    # original muscle align function
    #
    # start = time.time()
    # pc.bioMuscleAlign(in_dir, musl_path, outputF=out_dir)
    # end = time.time()
    # print(end - start)
    # # # ------------------------------
    # New multithreading based alignments
    #
    # start = time.time()
    # thread1 = pc.myThread(in_dir1, musl_path, outputF=out_dir)    
    # thread2 = pc.myThread(in_dir2, musl_path, outputF=out_dir)
    
    
    # thread1.start()
    # thread2.start()
    # end = time.time()
    # print(end - start)
    # for file_name in os.listdir(in_dir1):
    #     full_file_name = in_dir1 + file_name
    #     inputL.append(full_file_name)
    
    #print(inputL)
    #pc.bioMuscleAlignList(inputL, musl_path, out_dir)

    print(pc.gen_query_list(filePath))

    
if __name__ == '__main__':
    main()