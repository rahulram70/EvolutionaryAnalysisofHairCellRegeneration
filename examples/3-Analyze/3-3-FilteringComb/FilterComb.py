#!/usr/bin/env python3

import os
import sys



def main():
    # Creates the pathway to the ex_out_comb text file and the revised_ex_out_comb text file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
   
    # go up 2 directories
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))
    
    sys.path.append(res_dir)
    import procomp as pc

    start_dir = res_dir + "/resources/data-cleaned/combination_output/ex_out_comb.txt"
    revised_dir = res_dir + "/resources/data-cleaned/combination_output/revised_ex_out_comb.txt"
    removed_dir = res_dir + "/resources/data-cleaned/combination_output/removed_ex_out_comb.txt"
    
   

    # print(res_dir)
    # with open(in_fl, "r") as fl:
    #     fl = fl.read().splitlines()

    # Opens the ex_out_comb file and gives it reading permissions and opens the revised_dir file and gives it writing permissions
    ex_out_comb = open(start_dir, "r")
    revised_ex_out_comb = open(revised_dir , "w")
    removed_ex_out_comb = open(removed_dir, "w")
    pro_id = 0
    count = 0

    # Loops through all the lines in the ex_out_comb text file
    for line in ex_out_comb:
        line_list = []
    # Determines if the line is a zebrafish transcript header    
        if(line[0:7] == "ENSDART"):
            line_list = line.split()
    # Determines if there is either 0 non-regenerating or regenerating species and excludes it
            if(line_list[11] != "0" and line_list[13] != "0"):
    # Determines if there is a 1 to 1 comparison and excludes it            
                if(line_list[11] != "1" or line_list[13] != "1"):
    # Adds the header if it meets all above criteria. Also sets pro_id to 1 to allow the addition of its proteins to the new file                
                    revised_ex_out_comb.write(line)
                    pro_id = 1
                else:
                    removed_ex_out_comb.write(line)
                    pro_id = 0
                    count += 1
            else:
                removed_ex_out_comb.write(line)
                pro_id = 0
                count += 1
# Writes the proteins to the file if the transcript ID meets all criteria        
        elif(pro_id == 1):
            revised_ex_out_comb.write(line)
        else:
            removed_ex_out_comb.write(line)
                    
 

            
# Closes both files when done    
    ex_out_comb.close()
    revised_ex_out_comb.close()
    print(count)
    

if __name__ == '__main__':
    main()
