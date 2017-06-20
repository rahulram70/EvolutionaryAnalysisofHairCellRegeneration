#!/usr/bin/env python3

from __future__ import absolute_import
import os
import procomp as pc
import process as ps

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
    out_file = script_dir + "/ex_out_comb.txt"
    res_dir += "/resources/data-raw/transcript-ids/"
    out_file = open(out_file, "w+")

    L = []
    spid = [i.split()[0] for i in open(spid, "r").read().splitlines()]
    print(spid)

    for F in os.listdir(res_dir):
        if ("." not in F):
            ali_dir = res_dir + F + "/"
            L += pc.comb_TrPr(ali_dir, spid) 
    
    L = sorted(L)
    L = pc.comb_rm_dups(L, spid)
    for i in L:
        out_file.write(str(i) + "\n")
    #pc.comb_gen_combs(L, out_file)

    #for i in L:
    #    print (i)
    out_file.close()

if __name__ == '__main__':
    main()

    
