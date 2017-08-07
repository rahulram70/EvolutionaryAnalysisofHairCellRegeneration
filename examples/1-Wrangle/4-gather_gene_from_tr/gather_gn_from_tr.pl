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


#
# Set file paths to resources and output files
#
my $dir = getcwd();
my @dir_L = split "/", $dir;
$dir = join("/", @dir_L[0 ..$#dir_L - 3]);
$dir_res = "$dir/resources/data-raw/Transcript-ids-cleaned/";
$dir_out = "$dir/resources/data-raw/Gene-ids.txt";

opendir(DH, $dir_res);
my @dirs = readdir(DH);
close(DH);

# ------------------------
# Here we open the resource folder containing 
# all the species protein ids
#
foreach my $file (@dirs) {
    
    if (index($file, ".txt") != -1){
        print "Current file = $file\n";
        my $species = substr($file, 0, -4);

        # ---------------------------------
        # opens the file of ensembl protein ids
        #
        my $file_path = "$dir_res$file";

        
        open(my $fh, '<:encoding(UTF-8)', $file_path) or die "Could not open file '$file_path' $!";
        open(OUT, ">$dir_out");

        my $cur_tr = "";
        my $transcript_adaptor =
                    $registry->get_adaptor( 'Zebrafish', 'Core', 'Transcript' );
        my $cur_gn = "";
        while (my $line = <$fh>){
            my @line_L = split ' ', $line;
            if ( index($line_L[0], "ENS") != -1 && $line_L[0] ne $cur_tr ){
                $cur_tr = $line_L[0];
                #print "$cur_tr\n";
                my $transcript = $transcript_adaptor->fetch_by_stable_id($cur_tr);
                
                $cur_gn = $transcript->get_Gene->stable_id();
                printf OUT ("$cur_tr  $cur_gn\n");
            }
        }  
        close(OUT);
        close($fh);
        last;
    }
}
