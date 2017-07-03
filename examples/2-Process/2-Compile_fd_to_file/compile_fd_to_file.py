#!/usr/bin/env python3

from __future__ import absolute_import
import os
import procomp as pc
import process as ps

@pc.stats
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    inputFd = script_dir + "/testFd/"
    outputFile = script_dir + "/test.txt"
    
    ps.merge_folder_to_file(inputFd, outputFile)

if __name__ == '__main__':
    main()

    
