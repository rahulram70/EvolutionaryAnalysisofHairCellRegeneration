#!/usr/bin/env python2
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

def tr_to_pr(sp, tr):
    """
    OVERVIEW
        searches the database for a species (sp) and aquires
        the protein sequences associated with the provided
        transcript id (tr).
    INPUTS
        sp = species to search
        tr = transcript to convert to protein
    """

    return 0


print_species()

