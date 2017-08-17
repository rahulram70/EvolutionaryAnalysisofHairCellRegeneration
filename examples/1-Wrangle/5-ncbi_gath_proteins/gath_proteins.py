#!/usr/bin/env python3

import os
import sys
import time
from Bio import Entrez
from Bio import SeqIO

def main():
    """
    The goal of this script is to download all protein sequences for the 
    provided list of zebrafish genes.
    """

    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    # go up 2 directories to project dir
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))

    # Setup paths to needed directories
    #
    pl_fl = script_dir + "/gather_domain_info.pl"
    genes_path = res_dir + "/resources/data-raw/gene_names.txt"
    gene_output_path = res_dir + "/resources/data-raw-ncbi/Zebrafish_Protein_Sequences.txt"

    
    # Open Necessary files & setup variables
    #
    genes = open(genes_path, "r")
    genes_L = genes.read().splitlines()
    genes.close()
    gene_output = open(gene_output_path, "w+")  
    gene_output.write("Sequence Record ID, Gene, Protein Sequence\n")

    gene_cnt = len(genes_L[1:])
    cur_cnt = 1
    for i in range(1, len(genes_L) ):
        
        organism = "Danio rerio"
        stable_id = genes_L[i].split(", ")[1].split()[0]
        term = "{}[Gene] AND {}[Organism]".format(stable_id, organism)
        
        # ncbi login and query
        #   Always tell NCBI who you are
        Entrez.email = "harmonichemispheres@gmail.com"     

        # Search
        handle = Entrez.esearch(db="Protein", term=term)

        # Get search string to search for ids
        record = Entrez.read(handle)       
        ids = record['IdList']

        # setup entrez handle 
        handle2 = Entrez.efetch(db="Protein", id=ids, rettype="gb", retmode="text")

        # ------------------------------
        # Get zebrafish protein and id
        #     
        seq = ""
        longest_id = ""
        parse_handle = list(SeqIO.parse(handle2, "gb"))
        best_entry = ["", ""]

        for seq_record in parse_handle:        
            if "NP_" in seq_record.id:
                if len(seq_record.seq) > len(best_entry[1]) or "NP_" not in best_entry[0]:
                    best_entry[0] = seq_record.id
                    best_entry[1] = seq_record.seq
            
            # if no NP then look for XP
            #
            elif "NP_" not in best_entry[0]:
                if "XP_" in seq_record.id:
                    if len(seq_record.seq) > len(best_entry[1]) or "XP_" not in best_entry[0]:
                        best_entry[0] = seq_record.id
                        best_entry[1] = seq_record.seq
            
            # if no NP or XP then take longest other item
            #
            elif "NP_" not in best_entry[0] and "XP_" not in best_entry[0]:
                if "XP_" not in seq_record.id and "NP_" not in seq_record.id:
                    if len(seq_record.seq) > len(best_entry[1]):
                        best_entry[0] = seq_record.id
                        best_entry[1] = seq_record.seq
        
        print ("{}  BEST: {}   {}/{}".format(stable_id, best_entry[0], \
            cur_cnt, gene_cnt) )
        gene_output.write("{}, {}, {}\n".format(best_entry[0], stable_id, \
            best_entry[1]))
        
        # Close Entrez database handle
        handle.close()
        
        # keep track of query progress
        cur_cnt += 1
        
        # wait a few sec before next query to abide by ncbi policies
        time.sleep(2)

    gene_output.close()

if __name__ == '__main__':
    main()
