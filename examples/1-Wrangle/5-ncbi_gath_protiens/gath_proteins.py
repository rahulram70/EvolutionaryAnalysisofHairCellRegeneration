#!/usr/bin/env python3

import os
import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

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
    gene_out = res_dir + "/resources/data-raw/gene_list.txt"
    res_dir += "/resources/data-raw/protein-sequences/"


    #
    # Open Necessary files
    #
    genes = open(genes_path, "r")
    genes_L = genes.read().splitlines()

    #
    # Set up necessary variables
    #
    organism = "Danio rerio"
    stable_id = genes_L[3].split(", ")[1]
    term = "{}[Gene] AND {}[Organism]".format(stable_id, organism)
    print(term)

    #
    # ncbi login and query
    #
    Entrez.email = "harmonichemispheres@gmail.com"     # Always tell NCBI who you are

    # Search
    handle = Entrez.esearch(db="Protein", term=stable_id)   #, retmax='100'

    # Get search string to search for ids
    record = Entrez.read(handle)
    ids = record['IdList']
    print("{} hits".format(len(ids)))
    #for i in ids:
    #    print(i)

    # choose id to search
    seq_id = ids[0]

    # setup entrez handle 
    handle2 = Entrez.efetch(db="Protein", id=ids, rettype="gb", retmode="text")

    # ------------------------------
    # Get zebrafish protein and id
    #
    seq = ""
    for seq_record in SeqIO.parse(handle2, "gb"):
        if (organism in seq_record.annotations["source"]):
            print("{}_{}_len={}".format(seq_record.id, \
                seq_record.annotations["source"], len(seq_record)) )
            seq += "{}>{}".format(seq_record.annotations["source"], seq_record.seq)
            print(seq) 
            break
    
    # Close Entrez database handle
    handle.close()

if __name__ == '__main__':
    main()
