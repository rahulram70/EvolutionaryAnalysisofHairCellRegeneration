#!/usr/bin/env python3

import os
import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#
# Setup blast results object
#
class blast_res:
    def __init__(self, accession="", score=0, seq_len="", seq="", pos="", \
        quer="", subj="", e_val="", gaps="", org_query="", gene=""):
        self.accession = accession
        self.score = score
        self.seq_len = seq_len
        self.seq = seq
        self.pos = pos
        self.quer = quer
        self.org_query = org_query
        self.subj = subj
        self.e_val = e_val   
        self.gaps = gaps    
        self.gene = gene
    
    def __str__(self):
        ret_val = "ACC_ID: {}\nScore: {}\nSeq_Length: {}\nPositives: {}\nE_value: {}\nGaps: {}\nQuery: {}\n".format(self.accession, \
        self.score, self.seq_len, self.pos,self.e_val, self.gaps, self.quer)
        return ret_val

    def get_stats(self):
        ident = "ident: {}/{} ({}%)".format(len(self.seq), self.seq_len, \
            len(self.seq)/self.seq_len*100)
        positive = "Positives: {}/{} ({}%)".format(self.pos, self.seq_len, \
            self.pos/self.seq_len*100)
        gap = "Gaps: {}/{} ({}%)".format(self.gaps, self.seq_len, \
            self.gaps/self.seq_len*100)
        return {"ident": ident, "positive": positive, "gaps": gap}

    def query_cover(self):
        return len(self.subj)/len(self.quer)

    def stat(self):
        algn = len(self.seq)/self.seq_len*100
        ret_val = algn * self.org_query * self.query_cover()
        ret_val /= algn
        return ret_val
    
def parse_b_recs_algns(b_rec, thr):
    #
    # Parse out results
    #
    res_L = []
    for alignment in b_rec:
        for hsp in alignment.hsps:
            if hsp.expect < thr:
                
                # define variables from title line
                title_L = alignment.title.split("|")
                for i in title_L:
                    if "NP_" in i or "XP_" in i:
                        acc = i
                        break
                res = blast_res(accession=acc, score=int(hsp.score), seq_len=alignment.length, \
                  subj=hsp.sbjct, quer=hsp.query, e_val=hsp.expect, gaps=hsp.gaps, \
                  pos=hsp.positives)
                res_L.append(res)
    return res_L

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
    ncbi_xlsx = res_dir + "/resources/data-raw-ncbi/NCBI_Species_list.xlsx"
    res_dir += "/resources/data-raw/protein-sequences/"
    


    # 
    # Perform a Blast of the identified sequence
    #   
    #
    filter_str = "gallus gallus[ORGN]"
    org_seq = "NP_001315469.1"
    E_VALUE_THRESH = 0.04
    res_L = []
    
    #
    # FIRST Blast 
    #
    result_handle = NCBIWWW.qblast("blastp", "nr", org_seq, entrez_query=filter_str)
    blast_records = NCBIXML.read(result_handle)

    #
    # Parse out results
    #
    res_L = parse_b_recs_algns(blast_records.alignments, E_VALUE_THRESH)

    # Get protein with highest score
    #
    best = res_L[0]
    for i in res_L:
        print("acc=({}) score=({})".format(i.accession, i.score))
        if i.score > best.score:
            best = i
    print("BEST")
    print("acc=({}) score=({})".format(best.accession, best.score))
    
    #
    # Perform Reciprocal Blast
    #
    filter_str = "Danio rerio[ORGN]"
    result_handle = NCBIWWW.qblast("blastp", "nr", best.accession, entrez_query=filter_str)
    blast_records = NCBIXML.read(result_handle)
    
    #
    # Parse out results
    #
    del(res_L)
    res_L = parse_b_recs_algns(blast_records.alignments, E_VALUE_THRESH)
    
    # Check results for original Danio rerio protein
    #
    for i in res_L:
        print("acc=({}) score=({})".format(i.accession, i.score))
        if i.accession == org_seq:
            print("An ortholog has been found!!!")
            print(i)
    
if __name__ == '__main__':
    main()
