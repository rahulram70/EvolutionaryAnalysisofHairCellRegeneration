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
    align_path_dir = res_dir + "/resources/data-cleaned/aligned_seq/"
    spec_group_path = res_dir + "/resources/Species-Group.txt"

    # name, a_list, b_list = pc.comp_sort_region(align_path, spec_group_path)

    # print(name)
    # print(a_list)
    # print(b_list)
    print("Transcript ID")
    print("longest region")
    print("Largest region")
    for file in os.listdir(align_path_dir):
        
        align_path = align_path_dir + file
        name, a_list, b_list = pc.comp_sort_region(align_path, spec_group_path)
        print(name)
        print(a_list)
        print(b_list)


    
    # print(b_list)
    # print(c_list)



    
    #  print(hit)
    # print(count)
    # print(hit_list)
 

if __name__ == '__main__':
    main()