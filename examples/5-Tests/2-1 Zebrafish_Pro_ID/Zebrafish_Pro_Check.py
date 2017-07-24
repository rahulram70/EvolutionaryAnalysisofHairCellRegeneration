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
    temp_align_path = res_dir + "/resources/data-cleaned/transcript-alignments-all/"
    mstL_file_path = res_dir + "/resources/mstL_dictionary.txt"
    res_dir += "/resources/data-raw/Current-transcript-ids/" 
    L = []
    spid = [i.split() for i in open(spid_path, "r").read().splitlines()]
    print(spid)

    #  add all tr - pr associations to a list
    L = pc.comb_IdPr(tr_id_dir)
    
    #  sort tr - pr associations by tr id.
    L = sorted(L)

    #  remove any duplicates and generate table of all ortholog proteins
    L = pc.comb_rm_dups(L, list(zip(*spid))[0], ident="DART")

    #with open(script_dir + "/rm_dups.txt", "w+") as out:
    #    for i in L:
    #        out.write(str(i) + "\n")

    
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

    
    #print(mst_D)  
    mstL_file = open(mstL_file_path, "w+")    
    missing_pro_transIDs = []
    multiple_pro_transIDs = []
    for key in mst_D:
        count = 0
        for element in mst_D[key]:
            
            if("ENSDARP" in element):
                count += 1  
                   

                 
            
        if(count < 1 ):
            missing_pro_transIDs.append(key)
        if(count > 1):
            multiple_pro_transIDs.append(key)
        
    print(len(missing_pro_transIDs))
    print(multiple_pro_transIDs)
    for key, value in mst_D.items():
        mstL_file.write("{}\n".format(key))
        for line in value:
            mstL_file.write(str(line) + "\n")
    
    mstL_file.close()
    # mstL_file.write(i + "\n")
    #     for j in mstL[i]:
    #         for k in j:
    #             mstL_file.write(k + "\n")
 
if(__name__ == '__main__'):
    main()