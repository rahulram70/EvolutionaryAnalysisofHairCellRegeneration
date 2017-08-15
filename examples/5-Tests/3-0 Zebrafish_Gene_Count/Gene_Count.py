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

    gene_id_path = res_dir + "/resources/data-raw/Gene-ids.txt"
    gene_file = open(gene_id_path, "r")
    gene_count = 0
    gene_list = []

    for line in gene_file:
        if(line.split()[1][0:7] == "ENSDARG" and line.split()[1] not in gene_list):
            gene_count += 1
            gene_list.append(line.split()[1])

    print(gene_count)
    print(sorted(gene_list))





if(__name__ == '__main__'):
    main()