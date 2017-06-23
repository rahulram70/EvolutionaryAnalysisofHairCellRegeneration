#!/usr/bin/env python2
from __future__ import absolute_import





if __name__ == '__main__':
    import os
    import procomp as pc

    @pc.stats
    def main():
        
        script_dir = os.path.dirname(os.path.abspath(__file__))

        """ Neccesary to log into the ensembl database """
        import os
        from cogent.db.ensembl import HostAccount
        if 'ENSEMBL_ACCOUNT' in os.environ:
            host, username, password = os.environ['ENSEMBL_ACCOUNT'].split()
            account = HostAccount(host, username, password)
        else:
            account = None

        """ gathers the transcript id and protein sequence from gene """
        sp = "zebrafish"
        gn = "ENSDARG00000027279"
        from cogent.db.ensembl import Genome
        specie = Genome(Species=sp, Release="81", account=None)
        gene = specie.getGeneByStableId(StableId=gn)
        for tr in gene.Transcripts: 
            print(tr.StableId)
            print(tr.Cds)
    
    main()

    
