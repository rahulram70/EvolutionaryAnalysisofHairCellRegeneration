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
    mstL_file_path = res_dir + "/resources/mstL_dictionary.txt"
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
    #print(L)
    mstL_file = open(mstL_file_path, "w+")
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
    
    #print(mstL)      
    missing_pro_transIDs = []
    multiple_pro_transIDs = []
    for key in mstL:
        count = 0
        count_space = 0
        for value in mstL[key]:
            for element in value:
                if(element[0:7] == "ENSDARP"):
                    if(element[7:12] == "00000"):
                        count += 1  
                    elif(element[7:12] == "     "):
                        count_space += 1
                 
            
        if(count_space > 0 ):
            missing_pro_transIDs.append(key)
        if(count > 1):
            multiple_pro_transIDs.append(key)
        
    print(missing_pro_transIDs)
    print(len(multiple_pro_transIDs))
    for key, value in mstL.items():
        mstL_file.write("{}\n".format(key))
        for line in value:
            mstL_file.write(str(line) + "\n")
    
    mstL_file.close()
        ''' mstL_file.write(i + "\n")
        for j in mstL[i]:
            for k in j:
                mstL_file.write(k + "\n")
 '''
if(__name__ == '__main__'):
    main()