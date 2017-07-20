#!/usr/local/bin/perl
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',  # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

my $ID = 'ENSDARP00000131636';


my $seqmember_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','SeqMember');

# fetch a Member
my $seqmember = $seqmember_adaptor->fetch_by_stable_id($ID);
print $ID,"\n";
print $seqmember->sequence(),"\n";



#my $seq_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'SeqMember');
#my $seq_member = $seq_member_adaptor->fetch_by_stable_id('ENSDARP00000000005');

#my $family_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Family');
#my $protein_adaptor = $family_adaptor->fetch_by_SeqMember($seq_member);

#foreach $pr (@{$protein_adaptor->get_all_Members()}) {
#    my $id = $pr->stable_id();
#    if (index($id, "ENS") != -1) {
#        print $id . "\n";
#    }
    
#}
#print $protein_adaptor->stable_id() . "\n";
#print $protein_adaptor->description() . "\n";
#my $families = $family_adaptor->fetch_all_by_GeneMember($seq_member);


