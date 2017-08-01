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
    
    out_dir = script_dir + "/post-align/"
    in_dir1 = script_dir + "/pre_align1/"
    in_dir2 = script_dir + "/pre_align2/"
    in_dir = res_dir + "/resources/data-cleaned/transcript-longest-prealign/"
    musl_path = res_dir + "/resources/muscle/muscle3.8.31_i86win32.exe"
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

    threadL = pc.thread_creator(in_dir, 10, out_dir, musl_path)

    for obj in threadL:
        obj.start()
    
if __name__ == '__main__':
    main()