#!/usr/bin/env python3
from __future__ import absolute_import

import os

def main():
    """
    OVERVIEW:
        here we analyze transcript id resources gathered from ensembl
        with ortholog protein information to make sure the dataset is 
        consistent.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
    
    # go up 2 directories
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))


    spid = res_dir + "/resources/Species-Group.txt"
    spid = open(spid, "r").read().splitlines()
    trans_dir = res_dir + "/resources/data-raw/Current-transcript-ids/"
    master_tr_L = res_dir + "/resources/data-cleaned/Jiang_and_Steiner Merged Transcript list-out.csv"
    fl_out = script_dir + "/results.txt"

    # STEP 1 ---------------------------------
    # Setup a list of transcript ids from the master list
    #
    mstr_L = open(master_tr_L, "r").read().splitlines()
    mstr_L_tr = [i.split(',')[1] for i in mstr_L[1:]]
    

    # STEP 2 ---------------------------
    # Add all transcript for a species to a list
    # and count all files in the species folder
    #  
    for spe_name in os.listdir(trans_dir):
        if ("." not in spe_name):
            ali_dir = trans_dir + spe_name + "/"
            tr_L = []


            # STEP 3 -----------------------
            # generate a list of all unique tr ids for
            # this species
            #
            file_count = 0
            local_spe_L = []
            for file in os.listdir(ali_dir):
                
                if file.endswith(".txt"):    
                    fileloc = ali_dir + file
                    fileName = "" + file
                    fileString_split = open(fileloc, "r").read().splitlines()
                    file_count += 1
                    
                    for line in fileString_split[1:]:
                        if line.split()[0] not in local_spe_L: 
                            local_spe_L.append(line.split()[0])        
            
            
            # STEP 4 -----------------------
            # compare local list with master list to verify
            # all tr ids are present
            # 
            id_count = 0
            found = []

            for tr in mstr_L_tr:
                if tr in local_spe_L:
                    id_count += 1
                    found.append(tr)
                #
                # Use this to print out which ids are not present in the local
                # list, that ARE present in the master list.
                #
                #else:
                    #print("{} was not found in ({})".format(tr, spe_name))
                    #print(tr)
                    
            print("({}) files with ({}/{}) found ids in master for {}".format(file_count, id_count,len(mstr_L_tr), spe_name) )
            
            #
            # use this to see which ids in the local list
            # are not in the master list
            #
            #for i in local_spe_L:
            #    if i not in found:
            #        print("     {} was not found in local list".format(i))


if __name__ == '__main__':
    main()
