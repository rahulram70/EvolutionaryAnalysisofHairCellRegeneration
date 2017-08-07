#!/usr/local/bin/perl
use Bio::EnsEMBL::Registry;
use Cwd;
use File::Basename;
use Term::ProgressBar 2.00;


# --------------------------
# This sets up the ensembl registry and logs into the 
# ensembl mysql database.
#
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',  # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

# --------------------------------------------
# Set file paths to resources and output files
#
my $dir = getcwd();
my @dir_L = split "/", $dir;
$dir = join("/", @dir_L[0 ..$#dir_L - 3]);
$dir_res = "$dir/resources/data-raw/gene_list.txt";
$dir_out = "$dir/resources/data-raw/gene_names.txt";

# ------------------
# Open Needed Files
#
open(my $fh, '<:encoding(UTF-8)', $dir_res) or die "Could not open file '$dir_res' $!";
open(OUT, ">$dir_out");


# Setup ensembl database adaptor and variables
#
my $gene_adaptor =
    $registry->get_adaptor( 'Zebrafish', 'Core', 'Gene' );


#
# Iterate over genes
#
printf OUT ("ENS_stable_id, External_name\n");
while (my $row = <$fh>){
    chomp $row;
    my @line = split " ", $row;
    if (index($line[0], "G000") != -1){
        print "$line[0]\n";
        my $gene = $gene_adaptor->fetch_by_stable_id($line[0]);
        printf OUT ("%s, %s\n",
            $gene->stable_id(),
            $gene->external_name()
        );
    }
}
close($fh);
close(OUT);

