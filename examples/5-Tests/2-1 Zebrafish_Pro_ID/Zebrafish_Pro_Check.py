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

    # Setup paths to needed directories
    #
    spid_path = res_dir + "/resources/Species-Group.txt"
    out_file =  res_dir + "/resources/data-cleaned/"
    seq_dir = res_dir + "/resources/data-raw/protein-sequences/"
    temp_align_path = res_dir + "/resources/data-cleaned/transcript-alignments-all/"
    res_dir += "/resources/data-raw/Current-transcript-ids/"

    L = []
    spid = [i.split() for i in open(spid_path, "r").read().splitlines()]
    print(spid)

    #  add all tr - pr associations to a list
    for F in os.listdir(res_dir):
        if ("." not in F and F != "Icon\r"):
            ali_dir = res_dir + F + "/"
            L += pc.comb_IdPr(ali_dir) 
    
    #  sort tr - pr associations by tr id.
    L = sorted(L)

    #  remove any duplicates and generate table of all ortholog proteins
    L = pc.comb_rm_dups(L, list(zip(*spid))[0], ident="DART")
    

    #  Assemble 2d lists for each transcript
    reg = L[0][0]
    mstL = {}

    tmpL = []
    for i in L:
        if "DART" in i[0] and i[0] != reg :
            reg = i[0]
            mstL[ i[0] ] = list(tmpL)
            tmpL.clear()
        tmpL.append(i[1:])
        if "DART" in i[0] and L.index(i) == len(L)-1:
            reg = i[0]
            mstL[ i[0] ] = list(tmpL)
            tmpL.clear()
    
    
    #  Assemble dictionary of proteins ids and sequences
    pr_out_D = {}
    rng = 0
    for key, value in mstL.items():
        rng += 1
        pr_out_L = []
        for pr_L in value:
            for pr in pr_L:
                if (" " not in pr and "-" not in pr):
                    pr_out_L.append("{}".format(pr))
        pr_out_D[key] = pr_out_L
        if rng == 5:
            break
    