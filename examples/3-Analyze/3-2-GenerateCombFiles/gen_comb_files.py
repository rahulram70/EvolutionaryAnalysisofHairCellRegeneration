#!/usr/bin/env python3

import os
import procomp as pc

@pc.stats
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    in_fl = res_dir + "/3-0-CreatingCombinations/ex_out_comb.txt"
    # go up 2 directories
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))

    pr_seq_dir = res_dir + "/resources/data-raw/protein-sequences/"
    out_dir = res_dir + "/resources/data-cleaned/unaligned_seq/"

    log = pc.comb_gen_pro_files(in_fl, pr_seq_dir, out_dir)

    for i in log:
        print(i)

if __name__ == '__main__':
    main()
