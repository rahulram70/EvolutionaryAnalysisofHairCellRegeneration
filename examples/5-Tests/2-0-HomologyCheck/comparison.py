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
    align_path = res_dir + "/resources/data-cleaned/transcript-alignments-aligned/"
    test_file = res_dir + "/resources/data-cleaned/transcript-alignments-aligned/ENSDART00000000042.fasta"
    spid_path = res_dir + "/resources/Species-Group.txt"
    log_path = script_dir + "/log.csv"
    
    log = open(log_path, "w+")
    log.write("#,Transcript, R protein #, NR protein #, Total Proteins, \
    R >= 0.5 Sim, NR >= 0.5 sim, median sim, R species #, NR species #\n")
    c = 0
    # get average similarity for R and NR groups
    spid_tb = pc.spid_tb_gen(spid_path)

    for file in os.listdir(align_path):
        if ".fasta" in file:
            c += 1
            tr = file.split(".")[0] 
            tr_path = align_path + file
            print(tr)

            # generate log for percent similarity with zebrafish for transcript
            L = pc.comp_for_similarity(tr_path, "ENSDAR")


            g_r = 0
            spe_c_r = []
            spe_c_nr = []
            g_r_c = 0
            g_nr = 0
            g_nr_c = 0
            thr = 0.5
            median_L = [i[1] for i in L]
            median_L.sort()
            median = median_L[ len(median_L)//2 ]
            for spe in L:
                
                spe_gr = pc.spid_tb_get_group(spe[0], spid_tb)
                if spe_gr == "R":   
                    if spe[0][:6] not in spe_c_r:
                        spe_c_r.append(spe[0][:6]) 
                    if spe[1] >= thr:
                        g_r += 1
                    g_r_c += 1
                else:
                    if spe[0][:6] not in spe_c_nr:
                        spe_c_nr.append(spe[0][:6]) 
                    if spe[1] >= thr:
                        g_nr += 1
                    g_nr_c += 1
            log.write( "{},{},{},{},{},{},{},{},{},{}\n".format(c, tr, g_r_c, \
                g_nr_c, g_r_c + g_nr_c, g_r, g_nr, median, len(spe_c_r), len(spe_c_nr)) )


if __name__ == '__main__':
    main()