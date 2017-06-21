#!/usr/bin/env python2

from __future__ import absolute_import
import os
#import procomp as pc

"""
OVERVIEW
    functions in this module are focussed on query of ensembl
    using the pycogent module to access the ensembl mysql database.
"""

def ensembl_login():
    """ Neccesary to log into the ensembl database """
    import os
    from cogent.db.ensembl import HostAccount
    if 'ENSEMBL_ACCOUNT' in os.environ:
        host, username, password = os.environ['ENSEMBL_ACCOUNT'].split()
        account = HostAccount(host, username, password)
    else:
        account = None

def print_species():

    ensembl_login()

    from cogent.db.ensembl import Species
    print(Species)
    return 0

#@pc.stats
def tr_to_pr(sp, gn):
    """
    OVERVIEW
        searches the database for a species (sp) and aquires
        the protein sequences associated with the provided
        transcript id (tr).
    INPUTS
        sp = species to search
        tr = transcript to convert to protein
    NOTES
        sp, tr
    """

    ensembl_login()
    from cogent.db.ensembl import Genome
    specie = Genome(Species=sp, Release="81", account=None)
    gene = specie.getGeneByStableId(StableId=gn)
    for tr in gene.Transcripts: 
        print(tr.StableId)
        print(tr.Cds)
    """
    genes = specie.getGenesMatching(Symbol='brca2')
    for gene in genes:
        if gene.Symbol.lower() == 'brca2':
            break

    brca2 = gene # so we keep track of this reference for later on
    transcripts = brca2.Transcripts
    print (brca2.Symbol)
    """

    return 0


tr_to_pr('zebrafish','ENSDARG00000027279')

