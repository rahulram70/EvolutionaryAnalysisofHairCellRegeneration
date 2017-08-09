#!/usr/bin/env python3

import os
import sys
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
    gene_output_path = res_dir + "/resources/data-raw/Zebrafish_Protein_Sequences.txt"
    res_dir += "/resources/data-raw/protein-sequences/"
     


    #
    # Open Necessary files
    #
    genes = open(genes_path, "r")
    genes_L = genes.read().splitlines()
    gene_output = open(gene_output_path, "w+")
    
    genes_L.remove(genes_L[0])
    #
    # Set up necessary variables
    # >{}, {}, {}, {}\n".format(seq_record.id, " ".join(accession_list), stable_id, seq_record.seq)
    gene_output.write("Sequence Record ID, Gene, Protein Sequence\n")
    for i in range(len(genes_L)):
        
        organism = "Danio rerio"
        stable_id = genes_L[i].split(", ")[1].split()[0]
        term = "{}[Gene] AND {}[Organism]".format(stable_id, organism)
        #print(term)
        print(stable_id)
        #
        # ncbi login and query
        #
        Entrez.email = "harmonichemispheres@gmail.com"     # Always tell NCBI who you are

        # Search
        handle = Entrez.esearch(db="Protein", term=term)   #, retmax='100'

        # Get search string to search for ids
        record = Entrez.read(handle)
        
        ids = record['IdList']
        #print("{} hits".format(len(ids)))
        #for i in ids:
        #    print(i)

        # choose id to search
        
        
        #print(ids)

        # setup entrez handle 
        
        
        
        handle2 = Entrez.efetch(db="Protein", id=ids, rettype="gb", retmode="text")

        

            # ------------------------------
            # Get zebrafish protein and id
            #
            
        seq = ""
        count = 0
        
        longest_id = ""
        parse_handle = list(SeqIO.parse(handle2, "gb"))
        count = 0
        for seq_record in parse_handle:        
            
            if ( "NP_" in seq_record.id):
                print(seq_record.id + " NP confirmed")
                count += 1
                
                    #print(seq_record)
                    #print(seq_record.annotations["source"])
                    #print(stable_id)              
                    #accession_list = seq_record.annotations["accessions"]
                    #print(">{}, {}".format(seq_record.id, ", ".join(accession_list)))
                    #print("len={}".format(len(seq_record.seq)))
                gene_output.write("{},  {}, {}\n".format(seq_record.id, stable_id, seq_record.seq))
                    #print(">{}, {}, {}\n".format(seq_record.id, stable_id, seq_record.seq))
                

        
        if(count == 0):
            longest_rec = 0
            print("No NR found")       
            for element in parse_handle:
                print(len(element.seq))
                print(longest_rec)
                if("XP" in element.id and len(element.seq)> longest_rec):
                    longest_id = element.id
                    longest_rec = len(element.seq)
            
            gene_output.write("{},  {}, {}\n".format(element.id, stable_id, element.seq))
            print(stable_id + ": " + longest_id)
        
            if(longest_rec == 0):
                gene_output.write("{},  {}, {}\n".format("Error", stable_id, "Error"))
                     
            

                


                
                
        
        # Close Entrez database handle
        handle.close()

if __name__ == '__main__':
    main()
