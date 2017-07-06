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


    # STEP 1 ---------------------------
    # Add all transcript to a large list
    #
    spid = res_dir + "/resources/Species-Group.txt"
    spid = open(spid, "r").read().splitlines()
    trans_dir = res_dir + "/resources/data-raw/Current-transcript-ids/"
    fl_out = script_dir + "/results.txt"

    master_id_L = []
    master_result_L = []
    for F in os.listdir(trans_dir):
        if ("." not in F):
            ali_dir = trans_dir + F + "/"
            tr_L = []
            for file in os.listdir(ali_dir):
                if file.endswith(".txt"):
                    
                    fileloc = ali_dir + file
                    fileName = "" + file
                    fileString_split = open(fileloc, "r").read().splitlines()
                    count = 0
                    temp_id_L = []
                    for line in fileString_split:
                        if 0 < len(line.split()) <= 2:
                            count += 1
                            master_id_L.append( [ line.split()[0], file.split()[0], F] )
  
  
    

    # STEP 2 -----------------------------------------
    # Sort the large list created in the previous step
    #
    master_id_L.sort(key=lambda x: x[0])


    # STEP 3 --------------------------------------
    # iterate through list adding same tr ID lines to a list
    # and when complete iterate through that list checking for 
    # inconsistancy
    #
    check = master_id_L[0][0]
    spe_L = []
    out = open(fl_out, "w+")
    for i in master_id_L:
        
        spe_L.append(i[1][:7])

        if check != i[0] or i is master_id_L[-1]:
            
            missing = 0
            missing_L = []
            obtained = 0
            
            for s in spid:
                if s.split()[0][3:7] not in spe_L:
                    missing_L.append(s.split()[0][3:7])
                    missing += 1
                else:
                    obtained += 1
            
            if obtained != 61:
                out.write("{} has ids missing:   {}\n".format(check,missing_L) )
                print("{} has ids missing:   {}".format(check,missing_L) )
                #print("    missing ids ({})".format(missing) )
                #print("    present ids ({})".format(obtained) )
                
            check = i[0]
            spe_L = []
    out.close() 

if __name__ == '__main__':
    main()

