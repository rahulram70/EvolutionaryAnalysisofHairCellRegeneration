#!/usr/local/bin/perl
use Bio::EnsEMBL::Registry;
use Cwd;
use File::Basename;
use Term::ProgressBar 2.00;



$num_args = $#ARGV + 1;
my $select_spe = "";
if ($num_args != 1) {
    # use conditional 
    print "ERROR: Must provide species to use!\n";
    exit;
}
else{
    $select_spe = $ARGV[0];
}

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
$dir_res = "$dir/resources/data-raw/protein-sequences/";
$dir_out = "$dir/resources/data-raw/Domain-Features/";

opendir(DH, $dir_res);
my @dirs = readdir(DH);
close(DH);

# ------------------------
# Here we open the resource folder containing 
# all the species protein ids
#
foreach my $file (@dirs) {
    if (index($file, ".txt") != -1 && $file eq $select_spe){
        print "Current file = $file\n";
        my $species = substr($file, 0, -4);

        # ---------------------------------
        # opens the file of ensembl protein ids
        #
        my $file_path = "$dir_res$file";

        # get file length
        open(my $fh, '<:encoding(UTF-8)', $file_path) or die "Could not open file '$file_path' $!";
        $file_length++ while <$fh>;
        close $file_path;
        
        
        # ---------------------------
        # this opens a file for each species to write the retrieved 
        # data to.
        #
        my $spe_output = "$dir_out$file";
        open(FH, ">$spe_output");
        open(my $fh, '<:encoding(UTF-8)', $file_path) or die "Could not open file '$file_path' $!";

        # -------------------------------
        # after both input and output files have been opened we query
        # ensembl for each protein and write the results to the corresponding
        # output file.
        #
        my $line_count = 0;
        print "count=$line_count & length=$file_length\n";
        my $progress = Term::ProgressBar->new({name => 'Ensembl Query', count => $file_length, ETA => 'linear',} );
            $progress->max_update_rate(1);
            my $next_update = 0;
        print "\n";

        # Setup ensembl database adaptor
        #
        my $translation_adaptor =
            $registry->get_adaptor( $species, 'Core', 'Translation' );
        
        # Query for protein ids
        #
        while (my $row = <$fh>) {          
            chomp $row;
            $line_count = $line_count + 1;
            if (substr($row, 0, 1) eq ">" && index($row, "P000") != -1){
                $next_update = $progress->update($line_count) if($line_count > $next_update);
                my @temp_L = split ">", $row;
                my $stable_id = $temp_L[1];
                
                my $translation = $translation_adaptor->fetch_by_stable_id($stable_id);
                my $pfeatures = $translation->get_all_ProteinFeatures();
                while ( my $pfeature = shift @{$pfeatures} ) {
                    my $logic_name = $pfeature->analysis()->logic_name();
                    printf FH (
                        "%s, %d-%d, %s, %s, %s\n",
                        $stable_id,
                        $pfeature->start(), $pfeature->end(), $logic_name,
                        $pfeature->interpro_ac(),
                        $pfeature->idesc()
                    );
                }
            }
        }
        close(DH);
        close(FH);
        print "\n";
        exit;
    } 
}
