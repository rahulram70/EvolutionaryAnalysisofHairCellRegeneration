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
    tr_id_dir = res_dir + "/resources/data-raw/Transcript-ids-cleaned/"
    seq_dir = res_dir + "/resources/data-raw/protein-sequences/"
    temp_align_path = res_dir + "/resources/data-cleaned/transcript-longest-prealign/"
    pro_seq_path = res_dir + "/resources/data-raw/protein-sequences/"
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
    mst_edit_D = {}
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
    
    pro_seq_D = pc.gen_seq_hash_tb(seq_dir)
    count = 0
    for key in mst_D:
        spec_list = []
        pro_id_list = []
        for value in mst_D[key]:
            if(value != ""):
                if(value[0:7] in spec_list):
                    pro_index = spec_list.index(value[0:7])
                    pro_seq = pc.get_seq(pro_id_list[pro_index], pro_seq_D)
                    if(pro_seq != -1):
                        pro_seq_len = len(pro_seq)
                    new_pro_seq = pc.get_seq(value, pro_seq_D)
                    new_pro_len = len(pc.get_seq(value, pro_seq_D))
                    if(new_pro_len > pro_seq_len):
                        pro_id_list[pro_index] = value
                        count += 1
                else:
                    spec_list.append(value[0:7])
                    pro_id_list.append(value)
        
        mst_edit_D[key] = pro_id_list

    #print(mst_edit_D)

    
    # for key in mst_D:
    #     #print("{}\n{}".format(key, val))
    #     full_temp_align_path = temp_align_path + key + ".txt"
    #     print(full_temp_align_path)
    #     for val in mst_D[key]:
    #         pc.list_to_fasta(val, pro_seq_D, full_temp_align_path)
        
    # -------------------------
    #  Make files to be aligned
    #

    #  Generate Hashtable of sequences
    # seq_tb = pc.gen_seq_hash_tb(seq_dir)

    for key, value in mst_edit_D.items():
        out_path = temp_align_path + key + ".txt"
        pc.list_to_fasta(value, pro_seq_D, out_path)
        print("currently Making file for {}".format(key)) 
    
if __name__ == '__main__':
    main()
