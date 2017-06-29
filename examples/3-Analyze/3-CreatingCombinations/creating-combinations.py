#!/usr/bin/env python3

from __future__ import absolute_import
import os
import procomp.procomp as pc

"""
Here we generate combinations of transcripts that correlate to
Zebrafish genes and output a file of combinations in a 
semi-readable format.
"""

@pc.stats
def main():
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
    res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))
    spid = res_dir + "/resources/Species-Group.txt"
    print(script_dir)
    out_file = script_dir + "/ex_out_comb.txt"
    res_dir += "/resources/data-raw/transcript-ids/"

    L = []
    spid = [i.split()[0] for i in open(spid, "r").read().splitlines()]
    print(spid)

    for F in os.listdir(res_dir):
        if ("." not in F):
            ali_dir = res_dir + F + "/"
            L += pc.comb_TrPr(ali_dir, spid) 
    
    L = sorted(L)

    L = pc.comb_rm_dups(L, spid)
    
    pc.comb_gen_combs(L, out_file, 500)
    

if __name__ == '__main__':
    main()

    
