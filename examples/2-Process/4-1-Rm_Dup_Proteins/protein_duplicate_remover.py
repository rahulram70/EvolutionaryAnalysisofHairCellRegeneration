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
    species_pro_ID_database_list = []
    species_pro_seq_database_list = []

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
    
    pro_id, pro_seq = get_proID_database()
    # Creates the pathway to the procomp Directory    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    inputFile = os.path.abspath(os.path.join(script_dir, os.pardir))
    for i in range(2):
        inputFile = os.path.abspath(os.path.join(inputFile, os.pardir))
    
    # Adds the additional pathway to the protein-sequences
    inputFile = inputFile + "/resources/data-raw/protein-sequences/"
    
    # Intializes lists containing an empty list to accommodate the empty list in get_proID.
    #Lists will store the lists of animal proteins IDs and sequences.
    
    count = 0
    # For loop iterates through all the files in protein-sequences. 
    # Also intializes pro_id and pro_seq that will hold protein information for a species.
    for filename in os.listdir(inputFile):
        file_species_path = inputFile + filename
        an_new_file = open(file_species_path, "w")
        for element in range(len(pro_id[count])):
            if(element != 0):
                if(pro_id[count][element] != pro_id[count][element - 1]):
                    an_new_file.write(">" + pro_id[count][element] + "\n")
                    an_new_file.write(pro_seq[count][element] + "\n")
                    
                
                    
            else:
                an_new_file.write(">" + pro_id[count][element] + "\n")
                an_new_file.write(pro_seq[count][element] + "\n")
        
        count += 1    

if __name__ == '__main__':
    main()
 