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
        print(count)
    
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

    print(tmpList)
    for i in range(len(tmpList)):
        alist[left + i] = tmpList[i]
        blist[left + i] = tmpList2[i]

#alist = [3, 2 , 1]
#blist = ["tri", "cir", "sqr"]
#sort(alist, blist, 0, 2)
#print(alist, blist)

def main():

    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    inputFile = os.path.abspath(os.path.join(script_dir, os.pardir))
    for i in range(2):
        inputFile = os.path.abspath(os.path.join(inputFile, os.pardir))
    

    inputFile2 = inputFile + "\resources\data-raw\protein-sequences\\"

    #inputFile = inputFile + "\resources\data-raw\protein-sequences"

    

    print(inputFile2)

    for filename in os.listdir(inputFile):
    
        file_list = []
        pro_id = []
        pro_seq = []
        an_file = open(os.path.join(inputFile, filename), "r")
        file_list = an_file.readlines()
        count = 0
        
        for i in file_list:
            if count % 2 == 0:
                pro_id.append(file_list[count])
                count = count + 1
            else:
                pro_seq.append(file_list[count])                
                count = count + 1                        
                                     
           
    print(len(file_list))
    print(len(pro_id))
    print(len(pro_seq))
    print(pro_id[0], pro_id[2])
    print(pro_seq[1], pro_seq[3])
    

main()
