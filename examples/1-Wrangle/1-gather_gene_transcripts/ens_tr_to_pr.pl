#!/usr/bin/env perl
#
#
# a perl attempt gather transcript ids and their
# associated protein ids.

use Bio::EnsEMBL::Registry;
use Cwd;
use Cwd 'abs_path';
use strict;
use warnings;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);
my $transcript_adaptor =
    $registry->get_adaptor( 'Zebrafish', 'Core', 'Transcript' );

#@files = </Users/themusicman/Projects/Python/procomp/resources/data-raw/transcript-ids/*>;
my $some_dir = "/Users/themusicman/Projects/Python/procomp/resources/data-raw/transcript-ids";
opendir(my $dh, $some_dir) || die "Can't open $some_dir: $!";
while (readdir $dh) {
    
    my $some_dir_down_1 = "$some_dir/$_";

    if (index($some_dir_down_1, ".") == -1) {
        #print "$some_dir_down_1 does NOT contain .\n";
        opendir(my $dh_d1, $some_dir_down_1) || die "Can't open $some_dir_down_1: $!";
        while (readdir $dh_d1) {

            print "$some_dir_down_1/$_\n";
            my $filename = "$some_dir_down_1/$_";
            if (index($filename, "/.") == -1) {
                
                open(my $fh, '<:encoding(UTF-8)', $filename) || die "Could not open file '$filename' $!";
                
                while (my $row = <$fh>) {
                    chomp $row;
                    my @word_list = split(" ", $row);
                    if (index($word_list[0], "DART") != -1) {
                        print $word_list[0] . " " . $word_list[1] . "\n";
                        
                        my $stable_id = $word_list[0];
                        
                        my $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id);

                        print $transcript->stable_id(), "\n";
                        print $transcript->translation()->stable_id(), "\n";
                        print $transcript->translate()->seq(),         "\n";
                    }
                }
            } 
        }
        closedir $dh_d1;
    } 
    else {
        print "ERROR: $some_dir_down_1 does contain .\n";
    }
}
closedir $dh;

#foreach $folder (@files) {

    #my @words = split('/', $folder);
    #print $folder . "\n";
    #my $species_folder = $folder . @words[-1] . "/\n";
    #print $species_folder;
    #print @words[-1] . "\n";

#    @species_folders = <$folder>;
#    foreach $file (@species_folders){
#        print $file . "\n";
#    }

    #print $file . "\n";
#}

#my $stable_id = 'ENSDART00000000198';

#my $transcript_adaptor =
#  $registry->get_adaptor( 'Zebrafish', 'Core', 'Transcript' );
#my $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id);

#print $transcript->stable_id(), "\n";
#print $transcript->translation()->stable_id(), "\n";
#print $transcript->translate()->seq(),         "\n";

#print $transcript->translation()->transcript()->stable_id(), "\n";

