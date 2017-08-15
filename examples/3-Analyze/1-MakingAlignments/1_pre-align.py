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
    fil_path = script_dir + "/filtered/"
    tr_id_dir = res_dir + "/resources/data-raw/Transcript-ids-cleaned/"
    seq_dir = res_dir + "/resources/data-raw/protein-sequences/"
    temp_align_path = script_dir + "/unalignedSeqs/" #res_dir + "/resources/data-cleaned/transcript-prealign/"
    res_dir += "/resources/data-raw/Current-transcript-ids/"

    L = []
    spid = [i.split() for i in open(spid_path, "r").read().splitlines()]

    #  add all tr - pr associations to a list
    L = pc.comb_IdPr(tr_id_dir)
    
    #  sort tr - pr associations by tr id.
    L = sorted(L)

    #  remove any duplicates and generate table of all ortholog proteins
    L = pc.comb_rm_dups(L, list(zip(*spid))[0], ident="DART")
    
    #  Assemble 2d lists for each transcript
    reg = "ENSDART00000000005"  
    mst_D = {}
    tmpL = []
    length = len(L)
    for i in L:
        for pid in i:
            if "DART" in pid or (i == L[length-1] and i[-1] == pid): 
                if pid != reg:
                    mst_D[ reg ] = list(tmpL)
                    reg = pid
                    tmpL.clear()
            elif "-" not in pid and " " not in pid:
                tmpL.append(pid)
    
    #for key, val in mst_D.items():
    #    print("{}\n{}".format(key, val))
    # -------------------------
    #  Make files to be aligned
    #

    #  Generate Hashtable of sequences
    seq_tb = pc.gen_seq_hash_tb(seq_dir)
    thr = 0.5
    comparator = "DARP"
    for key, value in mst_D.items():
        out_path = temp_align_path + key + ".txt"
        pc.list_to_fasta(value, seq_tb, out_path)
        log = pc.comp_for_length(out_path, comparator, thr, \
          out_fl="{}{}.txt".format(fil_path, key))
        print(i)
            
        print("currently Making file for {}".format(key)) 
    
if __name__ == '__main__':
    main()
