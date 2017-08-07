#!/usr/bin/env python3

import os
import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def main():
    """
    This script blasts each protien sequence from those downloaded in 
    gath_proteins.py and performs a reciprocal blast to identify the species
    that have the highest similarity to the zebrafish and there by ensure we
    are working with orthologs and simultaneously eliminate the possibility of
    paralogs.
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
    # Setup blast results object
    #
    class blast_res:
        def __init__(self, accession, score, seq_len, seq, pos, quer, e_val, gaps):
            self.accession = accession
            self.score = score
            self.seq_len = seq_len
            self.seq = seq
            self.pos = pos
            self.quer = quer
            self.e_val = e_val   
            self.gaps = gaps    
        
        def __str__(self):
            ret_val = "{}\n{}\n{}\n{}\n{}\n{}\n{}\n".format(self.accession, \
            self.score, self.seq_len, self.pos, self.quer, self.e_val, self.gaps)
            return ret_val

        def get_stats(self):
            ident = "ident: {}/{} ({}%)".format(len(self.seq), self.seq_len, \
                len(self.seq)/self.seq_len*100)
            positive = "Positives: {}/{} ({}%)".format(self.pos, self.seq_len, \
                self.pos/self.seq_len*100)
            gap = "Gaps: {}/{} ({}%)".format(self.gaps, self.seq_len, \
                self.gaps/self.seq_len*100)
            return {"ident": ident, "positive": positive, "gaps": gap}


    # 
    # Perform a Blast of the identified sequence
    #   
    #
    filter_str = "gallus gallus[ORGN]"
    seq = "NP_001315469.1"
    test = blast_res("new guy", 3, 5, 5, "j", 6, 8)
    print(test)
    """
    E_VALUE_THRESH = 0.04
    
    result_handle = NCBIWWW.qblast("blastp", "nr", seq, entrez_query=filter_str)
    blast_records = NCBIXML.read(result_handle)


    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('Score:', hsp.score)
                print('ident:', hsp.identities)
                print('query:', hsp.query)
                print('length:', alignment.length)
                print('e value:', hsp.expect, "\n")
                
                #print(hsp.query[0:75] + '...')
                #print(hsp.match[0:75] + '...')
                #print(hsp.sbjct[0:75] + '...')
    

    """

if __name__ == '__main__':
    main()
