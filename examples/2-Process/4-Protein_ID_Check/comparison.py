#!/usr/bin/env python3
import os


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

# alist = ["tri", "cir" , "sqr"]
# blist = ["tri", "cir", "sqr"]
# sort(alist, blist, 0, 2)
# print(alist, blist) 

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
# alist = ["tri", "cir" , "sqr"]
# blist = ["tri", "cir", "sqr"]
# sort(alist, blist, 0, 2)
# print(alist, blist)
# print(search(alist, "tri")) 

def get_proID():
    
    # Creates the pathway to the procomp Directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    inputFile = os.path.abspath(os.path.join(script_dir, os.pardir))
    for i in range(2):
        inputFile = os.path.abspath(os.path.join(inputFile, os.pardir))

    # Adds the additional pathway to the Current-transcript-ids
    inputFile = inputFile + "/resources/data-raw/Current-transcript-ids/"
    species_pro_ID_list = []

    # First for loop iterates through all the folders in Current-transcript-ids. 
    # Also creates an empty list that will store all the protein IDs for a species
    for subdir, dirs, files in os.walk(inputFile):              
        species_pro_ID = []
        # Second for loop iterates thorugh each textfile in each animal folder and opens it.
        for file in files:
            species_file = open(os.path.join(subdir, file))
            # Third for loop iterates through each line in the textfile.
            # Also, determines if a line contains a protein ID and adds the ID to the species_pro_ID list if it does.
            for line in species_file:
                if(len(line) == 38):
                    species_pro_ID.append(line[19:37])
        
        
        species_pro_ID_list.append(sorted(species_pro_ID))
    return species_pro_ID_list
    
   # print(species_pro_ID_list[1])

def get_proID_database():

    # Creates the pathway to the procomp Directory    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    inputFile = os.path.abspath(os.path.join(script_dir, os.pardir))
    for i in range(2):
        inputFile = os.path.abspath(os.path.join(inputFile, os.pardir))
    
    # Adds the additional pathway to the protein-sequences
    inputFile = inputFile + "/resources/data-raw/protein-sequences/"
    
    # Intializes lists containing an empty list to accommodate the empty list in get_proID.
    #Lists will store the lists of animal proteins IDs and sequences.
    species_pro_ID_database_list = [[]]
    species_pro_seq_database_list = [[]]

    # For loop iterates through all the files in protein-sequences. 
    # Also intializes pro_id and pro_seq that will hold protein information for a species.
    for filename in os.listdir(inputFile):
        pro_id = []
        pro_seq = []
        file_species_path = inputFile + filename
        an_file = open(file_species_path, "r")
        ID_Check = 0
        # Second for loop iterates through every line in a file and the if statements determine if the line is an ID or sequence.
        # Accommodates for various errors such as multiple sequences or IDs in a row.
        for line in an_file:
            if(line[0:4] == ">ENS"):
                if(ID_Check == 0):
                    pro_id.append(line[1:-1])
                    ID_Check = 1
                else:
                    pro_seq.append("")
                    pro_id.append(line[1:-1])
            elif(ID_Check == 1):
                pro_seq.append(line[1:-1])
                ID_Check = 0
            else:
                pro_id.append("")
                pro_seq.append(line[1:-1])
        # Sorts the pro_id list and pro_seq using the sort function.
        # Adds the sorted lists into the previously intialized lists and returns them.
        sort(pro_id, pro_seq, 0, len(pro_id) - 1 )
        species_pro_ID_database_list.append(pro_id)
        species_pro_seq_database_list.append(pro_seq)

    return species_pro_ID_database_list, species_pro_seq_database_list
    

   
def main():

    pro_ID = get_proID()
    pro_ID_database, pro_seq_database = get_proID_database()
    
   # print(pro_ID[0], pro_ID_database[0], pro_seq_database[0])
    
    
    error_proID_List = []
    error_proseq_List = []
    for List in range(len(pro_ID)):
        for element in pro_ID[List]:
            
            #print(pro_ID_database[List])
            #print(element)
            found_index = search(pro_ID_database[List], element)
            
            if (found_index != -1):
                if(pro_seq_database[List][found_index][-1:] != "*"):
                    error_proseq_List.append(element)
            else:
                error_proID_List.append(element)

    print(error_proseq_List)
    print(error_proID_List)
        
if __name__ == '__main__':
    main()
