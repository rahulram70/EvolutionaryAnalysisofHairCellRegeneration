#!/usr/bin/env python3

from __future__ import absolute_import
import os

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    inputFile = os.path.abspath(os.path.join(script_dir, os.pardir))
    inputFile = os.path.abspath(os.path.join(inputFile, os.pardir))

    inputFile = inputFile + "/resources/Species-Group.txt"

    
    in_list = open(inputFile, "r").read().splitlines()
    out_list = open(inputFile, "w")
    
    for i in sorted(in_list):
        out_list.write(i + "\n")

if __name__ == '__main__':
    main()
