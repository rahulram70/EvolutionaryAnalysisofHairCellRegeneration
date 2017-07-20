#!/usr/bin/env python3

import os



def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
   
    # go up 2 directories
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))
    
    res_dir = res_dir + "/resources/resources/data-cleaned/combination_output/ex_out_comb.txt"
    revised_dir = res_dir + "/resources/resources/data-cleaned/combination_output/revised_ex_out_comb.txt"
    
   

    print(res_dir)
    # with open(in_fl, "r") as fl:
    #     fl = fl.read().splitlines()

    ex_out_comb = open(res_dir, "r")
    revised_ex_out_comb = open(revised_dir , "w")
    

    for line in ex_out_comb:
        line_list = []
        
        if(line[0:3] == "ENS"):
            line_list = line.split()
            if(line_list[11] != 0 and line_list[13] != 0):
                if(!(line_list[11] == 1 and line_list[13] == 1)):
                   revised_ex_out_comb.write(line)
                   pro_id = 1
                else:
                    pro_id = 0
            else:
                pro_id = 0
        
        elif(pro_id == 1):
            revised_ex_out_comb.write(line)
                   


            
    
    ex_out_comb.close()
    
    

if __name__ == '__main__':
    main()
