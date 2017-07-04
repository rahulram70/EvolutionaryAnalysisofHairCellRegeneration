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

    # go up 2 directories
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))

    spid = res_dir + "/resources/Species-Group.txt"
    print(script_dir)
    out_file = script_dir + "/ex_out_comb.txt"
    res_dir += "/resources/data-raw/Protein-ids from transcript query/"

    L = []
    spid = [i.split()[0] for i in open(spid, "r").read().splitlines()]
    print(spid)

    # add all tr - pr associations to a list
    for F in os.listdir(res_dir):
        if ("." not in F):
            ali_dir = res_dir + F + "/"
            L += pc.comb_IdPr(ali_dir, spid) 
    
    # sort tr - pr associations by tr id.
    L = sorted(L)

    # remove any duplicates and generate table of all ortholog proteins
    L = pc.comb_rm_dups(L, spid, ident="DART")
    
    # generate combinations of orthologs with protein ids
    # transcripts with too many combinations will be refactored
    # to a more managable count.
    pc.comb_gen_combs(L, out_file, 500)
    

if __name__ == '__main__':
    main()

    
