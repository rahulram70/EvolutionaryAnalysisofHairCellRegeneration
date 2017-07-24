#!/usr/bin/env python3

import os
import sys

def sort(alist, blist,  left,  right):

    if right == left: return
    middle = (left + right) // 2
    sort(alist, blist, left, middle)
    sort(alist, blist, middle + 1, right)
    merge(alist, blist, left, middle, right)

def merge(alist, blist, left, middle, right):
    tmpList = []
    tmpList2 = []
    count = 0
    for i in range(right - left + 1):
        count += 1
        tmpList.append("")
        tmpList2.append("")
        
    
    index1 = left
    index2 = middle + 1
    indx = 0
    while(index1 <= middle and index2 <= right):
        if alist[index1] < alist[index2]:
            tmpList[indx] = alist[index1]
            tmpList2[indx] = blist[index1]
            index1 += 1
        else:
            tmpList[indx] = alist[index2]
            tmpList2[indx] = blist[index2]
            index2 += 1
        indx += 1

    while index1 <= middle:
        tmpList[indx] = alist[index1]
        tmpList2[indx] = blist[index1]
        index1+= 1
        indx+= 1

    while index2 <= right:
        tmpList[indx] = alist[index2]
        tmpList2[indx] = blist[index2]
        index2 += 1
        indx += 1
  
    for i in range(len(tmpList)):
        alist[left + i] = tmpList[i]
        blist[left + i] = tmpList2[i]
def search(alist, srchVal):
    lb = 0
    ub = len(alist) - 1
    while(lb <= ub):
        mid = (lb + ub) // 2
        if(alist[mid] == srchVal):
            return mid
        elif(srchVal > alist[mid]):
            lb = mid + 1
        else:
            ub = mid - 1
    return -1

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
    align_path = res_dir + "/resources/data-cleaned/transcript-prealign/"
    spid_path = res_dir + "/resources/Species-Group.txt"
    log_path = script_dir + "/log.csv"
    
    log = open(log_path, "w+")
    log.write("Transcript ID, # R Species, # NR Species, Below 50% in zebrafish length, Above 50% in zebrafish length, error\n")

    c = 0
    spid_tb = pc.spid_tb_gen(spid_path)
    for file in os.listdir(align_path):
        
        if ".fasta" in file:
            c += 1
            tr = file.split(".")[0] 
            tr_path = align_path + file
            #print(tr)

            # generate log for percent similarity with zebrafish for transcript
            #L = pc.comp_for_similarity(tr_path, "ENSDAR")

            # get average similarity for R and NR groups
            
            spec_pro_id_L = []
            
            pro_seq_L = []
            pro_id_L = []
            regen_spec = 0
            non_regen_spc = 0
            upper = 0
            lower = 0
            error = ""
            
            #median_L = [i[1] for i in L]
            #median_L.sort()
            #median = median_L[ len(median_L)//2 ]
            file_species_path = align_path + file
            an_file = open(file_species_path, "r")
        
            for line in an_file:
                if(line[0:4] == ">ENS"):
                    
                    pro_id_L.append(line[1:8])
                    spec_pro_id_L.append(line)
                    
                    
                else:
                    pro_seq_L.append(line)
                    spec_pro_id_L.append(line)
                
            if(len (pro_id_L) > 1):
                sort(pro_id_L, pro_seq_L, 0, len(pro_id_L) - 1)
            
            zebrafish_index = search(pro_id_L, "ENSDARP")
            
            if(zebrafish_index != -1):
                error = ""
                zebrafish_len = len(pro_seq_L[zebrafish_index])
                for line in spec_pro_id_L:
                    
                    if(line[0:4] == ">ENS"):
                        pro_id = line[1:]
                        #print(pc.spid_tb_get_group(pro_id, spid_tb))
                        if( pc.spid_tb_get_group(pro_id, spid_tb) ==  "R"):
                            regen_spec += 1
                        elif(pc.spid_tb_get_group(pro_id, spid_tb) == "nr" or "NR"):
                            non_regen_spc += 1

                    elif(len(line) >= zebrafish_len):
                        upper += 1
                    else:
                        lower += 1
            else:
                error = "Does not contain zebrafish protein"
                for line in spec_pro_id_L:
                    if(line[0:4] == ">ENS"):
                                       
                        if( pc.spid_tb_get_group(pro_id, spid_tb) == "r" or "R"):
                            regen_spec += 1
                        elif(pc.spid_tb_get_group(pro_id, spid_tb) == "nr" or "NR"):
                            non_regen_spc += 1

            log.write( "{},{},{},{},{},{}\n".format(file[:-6], " " + str(regen_spec), " " + str(non_regen_spc), " " + str(lower), " " +str(upper),  " " + error))
            #break            
            
            


            

if __name__ == '__main__': 
    main()
