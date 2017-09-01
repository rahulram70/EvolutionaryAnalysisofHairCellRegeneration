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
    def __init__(self, specie="", accession="", score=0, seq_len="", pos="", \
        quer="", subj="", e_val="", gaps="", org_query="", gene=""):
        self.accession = accession
        self.specie = specie
        self.score = score
        self.seq_len = seq_len
        self.pos = pos
        self.quer = quer
        self.org_query = org_query
        self.subj = subj
        self.e_val = e_val   
        self.gaps = gaps    
        self.gene = gene
    
    def __str__(self):
        ret_val = "Species: {},  ACC_ID: {},  Score: {},  E_val: {}".format(self.specie, \
              self.accession, self.score, self.e_val)
        return ret_val

    def get_stats(self):
        ident = "ident: {}/{} ({}%)".format(len(self.seq), self.seq_len, \
            len(self.seq)/self.seq_len*100)
        positive = "Positives: {}/{} ({}%)".format(self.pos, self.seq_len, \
            self.pos/self.seq_len*100)
        gap = "Gaps: {}/{} ({}%)".format(self.gaps, self.seq_len, \
            self.gaps/self.seq_len*100)
        return {"ident": ident, "positive": positive, "gaps": gap}
    
def parse_b_recs_algns(b_rec):
    
    # Parse out results
    #
    ret_val = {}
    for alignment in b_rec:
        for hsp in alignment.hsps:

            sp = alignment.title.split("[")[1].split("]")[0]

            print("    title:::" + alignment.hit_def)
            # define variables from title line
            acc = alignment.hit_id.split("|")[-2]
            res = blast_res(accession=acc, score=int(hsp.score), \
                subj=hsp.sbjct, e_val=hsp.expect, specie=sp)
            
            if sp in ret_val:
                ret_val[sp].append(res)
            else:
                ret_val[sp] = []
                ret_val[sp].append(res)

    return ret_val


def main():
    """
    This script blasts each protien sequence from those downloaded in 
    gath_proteins.py and performs a reciprocal blast to identify the species
    that have the highest similarity to the zebrafish and there by ensure we
    are working with orthologs and simultaneously eliminate the possibility of
    paralogs.
    """

    E_VALUE_THRESH = 0.001
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
    dan_rer_proteins = res_dir + "/resources/data-raw-ncbi/danio_rerio_proteins.txt"
    
    # 
    # Perform a Blast of the identified sequence
    #   
    sp = '"Macaca mulatta"[Organism]' # OR Cercocebus atys[Organism]
    dan_seqs = ("MRARYMLGFGVLCLLTWAVYVSESSPTDKQAKLRLMSVRRARQNKSRQGQIIKTRASSTWINDKRQQMQGDEGVQFRRNVQSRKVLAQRRPAQGGTRNPIQQDDGTPGARARVSRMPSSAGSPNLLASFAGKNRLLVISAPHDSDGYYRLMMSLLKPEVYCELAERHVHQIVMFHQEGELGGKIRRITNEGKVMEEPLDTALIPRLMTFLKLEKGKFGMVLLKKTLQVEERYPYPVRLEAMYEVIDQSPMRKMEKVRQKGFVQKCKGAGVEGQVVEGVLTTDSQVDPPTERRKEIRKPIRRPTTTTTPAPTRPTTTTTTTKATTTTTTRPPTTTRSTTTTTTTTTTTTRPTTTTTRTTTTPRTTRANTTPQWIPAHKTTAEPYYYNRRDRYQTTSPPTDSARYRDNHTSKKEYNHRHTNTIPTQHKPTKVRPTKKKNGDKDISNAYEEKYDVGVPTDAYPEEKEEEIVPTKRGKGKTDKKKKKDKTDKLSKKDKAERRGKDGKGGKKNGKKVPKILEKEDYQKPTKRPPPPPPPKGTLATFLDYFESRRRLILITSPTEENSMYIQQRDEYLEHVCEMAIRKVTIITIFGTFRNSTMKIDHYQLEKDKPMKGLRQEDLENQDLIMELRKEYGMTYNDFYVVLTDLDMKAKQYYEVPIAMKAVFDYIDTFSSRIREMEQQKRDGVTCKKEDKPRSLENFLSRFRWRRRLFVISAPNDEEWAYQQQLYALTSQACNLGLRHVSVLKLVGTDLLDMGGVLELYPINGSATVEREGISATLVRDIRNYFQISPEYFSMLLVGKDGNVKSWYPSPMWSMAIIYDLIDSMQLRRQEMAIQQSLGMRCPEDEYGGYGYHHHEGYQEGYHQGYGY", \
                "NP_001037778.2") 
    res_L = []
    
    for org_seq in dan_seqs:

        #
        # FIRST Blast 
        #
        print("-- QUERY: {}, {} --".format(org_seq, sp))
        result_handle = NCBIWWW.qblast("blastp", "nr", org_seq, \
          entrez_query=sp, expect=E_VALUE_THRESH, \
          ncbi_gi=True)
          #, hitlist_size=5
        print("S2")
        blast_records = NCBIXML.read(result_handle)

        #
        # Parse out results
        #
        res_d = parse_b_recs_algns(blast_records.alignments)

        # Get protein with highest score
        #
        if len(res_d) > 0:
            best_L = []
            for key, val in res_d.items():
                best = blast_res(score=0)
                for hit in val:
                    if hit.score > best.score:
                        best = hit
                best_L.append(best)
                print("BEST || " + str(best))

            true_orthologs_L = []
            for sp in best_L:
                
                # Perform Reciprocal Blast
                #
                result_handle = NCBIWWW.qblast("blastp", "nr", best.accession, \
                  entrez_query="Danio rerio[organism]", hitlist_size=5, \
                  expect=E_VALUE_THRESH) 
                blast_records = NCBIXML.read(result_handle)

                # Parse out results
                #
                del(res_d)
                res_d = parse_b_recs_algns(blast_records.alignments)

                # Check results for original Danio rerio protein
                #
                for i in res_d["Danio rerio"]:
                    if i.accession == org_seq:
                        print("An ortholog has been found!!!")
                        print("   RESULT: " + str(i))
                        true_orthologs_L.append(sp)
                        break
            
            for i in true_orthologs_L:
                print("SHOW: " + str(i))
        print("\n\n")
    
if __name__ == '__main__':
    main()
