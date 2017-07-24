#!/usr/bin/env python3

import os
import sys

def comb_all_IdPr(id_dir, name_filter=""):
    """
    OVERVIEW:
        combines a folder of Id -- protein associations
        into a long list. Id are typically (XXXG, XXXT, XXXP)
    INPUTS:
        alignments_dir = path to folder with alignments
        spgr           = string for species id to look for (i.e. "DARP")
        name_filter    = you can choose to use only files with the given
                         string filter, default set to accept all files.
        print_computed = will print out the files used.
    """
    tr_pr_L = []
    for file in os.listdir(id_dir):
        if file.endswith(".txt") and name_filter in file:
            fileloc = id_dir + file
            fileName = "" + file
            fileString = open(fileloc, "r")
            fileString_split = fileString.read().splitlines()
            
            for line in fileString_split:
                tr_pr_L.append(str(line))
            fileString.close()
    print("comb_TrPr() has finished")
    return tr_pr_L

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
    out_dir =  res_dir + "/resources/data-raw/Transcript-ids-cleaned/"
    temp_align_path = res_dir + "/resources/data-cleaned/transcript-alignments-all/"
    res_dir += "/resources/data-raw/Current-transcript-ids/"

    L = []

    #  add all tr - pr associations to a list
    for F in os.listdir(res_dir):
        if ("." not in F and F != "Icon\r"):
            ali_dir = res_dir + F + "/"
            
            out = out_dir + F + ".txt"
            out_fl = open(out, "w+")

            L = comb_all_IdPr(ali_dir) 
            out_fl.write(L[0] + "\n")

            for line in L:
                ln_L = line.split()
                if len(ln_L) == 2 and "T00" in ln_L[1]:
                    pass
                elif "stable ID" not in line:
                    out_fl.write( "{}\n".format(line) )
            out_fl.close()
            


if __name__ == '__main__':
    main()