#!/usr/bin/env python3

from __future__ import absolute_import
import os
import procomp.procomp as pc
import statistics as sts

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
            L += pc.comb_IdPr(ali_dir, spid, name_filter=" 1-500") 

    L = sorted(L)
    L = pc.comb_rm_dups(L, spid, ident="DART")
    L_Tr = pc.comb_gen_combs(L, out_file, 500, ident="DART", refac=0, w=0)
    
    
    #
    # Get Protein results
    #    
    print("======== PROTEIN ========")
    L = []
    for F in os.listdir(prot_dir):
        if ("." not in F):
            ali_dir = prot_dir + F + "/"
            L += pc.comb_IdPr(ali_dir, spid) 
 
    L = sorted(L)
    L = pc.comb_rm_dups(L, spid, ident="DARP") 
    L_Pr = pc.comb_gen_combs(L, out_file, 500, ident="DARP", refac=0, w=0)
    
    
    #
    # Get Gene results
    #    
    print("======== GENE ========")
    L = []
    for F in os.listdir(gene_dir):
        if ("." not in F):
            ali_dir = gene_dir + F + "/"
            L += pc.comb_IdPr(ali_dir, spid) 
 
    L = sorted(L)
    L = pc.comb_rm_dups(L, spid, ident="DARG")
    L_Gn = pc.comb_gen_combs(L, out_file, 500, ident="DARG", refac=0, w=0)

    
    #
    # Visualize Data
    #
    print("\n\n======== RESULTS ========")
    L_Gn_comb = [int(i.split()[-1]) for i in L_Gn]
    L_Tr_comb = [int(i.split()[-1]) for i in L_Tr]
    L_Pr_comb = [int(i.split()[-1]) for i in L_Pr]
    print("Gene        Average Combinations:", sum(L_Gn_comb)/len(L_Gn_comb) )
    print("Transcript  Average Combinations:", sum(L_Tr_comb)/len(L_Tr_comb) )
    print("Protein     Average Combinations:", sum(L_Pr_comb)/len(L_Pr_comb) )
    
    print("\nGene        Median Combinations:", sts.median(L_Gn_comb) )
    print("Transcript  Median Combinations:", sts.median(L_Tr_comb) )
    print("Protein     Median Combinations:", sts.median(L_Pr_comb) )

    print("======== ======== ========\n\n")

    
if __name__ == '__main__':
    main()

