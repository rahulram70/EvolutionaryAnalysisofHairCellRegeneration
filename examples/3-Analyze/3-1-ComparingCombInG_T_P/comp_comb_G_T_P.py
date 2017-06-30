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
    out_file = script_dir + "/ex_out_comb.txt"
    trans_dir = res_dir + "/resources/data-raw/Protein-ids from transcript query/"
    prot_dir = res_dir + "/resources/data-raw/Protein-ids from protein query/"
    gene_dir = res_dir + "/resources/data-raw/Protein-ids from gene query/"

    spid = [i.split()[0] for i in open(spid, "r").read().splitlines()]
    print(spid)
    
    #
    # Get transcript results
    #
    
    print("======== TRANSCRIPT ========")
    L = []
    for F in os.listdir(trans_dir):
        if ("." not in F):
            ali_dir = trans_dir + F + "/"
            L += pc.comb_TrPr(ali_dir, spid, name_filter=" 1-500") 

    L = sorted(L)
    for i in L:
        print(i)
    L = pc.comb_rm_dups(L, spid, ident="DART")
    pc.comb_gen_combs(L, out_file, 500, ident="DART", refac=0)
    
    
    #
    # Get Protein results
    #    
    """
    print("======== PROTEIN ========")
    L = []
    for F in os.listdir(prot_dir):
        if ("." not in F):
            ali_dir = prot_dir + F + "/"
            L += pc.comb_TrPr(ali_dir, spid) 
            #print (L)
 
    L = sorted(L)
    for i in L:
        print(i)
    L = pc.comb_rm_dups(L, spid, ident="DARP") 
    pc.comb_gen_combs(L, out_file, 500, ident="DARP", refac=0)
    """
    
    #
    # Get Gene results
    #    
    """
    print("======== GENE ========")
    L = []
    for F in os.listdir(gene_dir):
        if ("." not in F):
            ali_dir = gene_dir + F + "/"
            L += pc.comb_TrPr(ali_dir, spid) 
 
    L = sorted(L)
    
    L = pc.comb_rm_dups(L, spid, ident="DARG")
    
    pc.comb_gen_combs(L, out_file, 500, ident="DARG", refac=0)
    """
if __name__ == '__main__':
    main()

