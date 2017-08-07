#!/usr/local/bin/perl
use Bio::EnsEMBL::Registry;
use Cwd;

# --------------------------
# This sets up the ensembl registry and logs into the 
# ensembl mysql database.
#
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',  # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);
my $seqmember_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','SeqMember');

my $dir = getcwd();
my @dir_L = split "/", $dir;
$dir = join("/", @dir_L[0 ..$#dir_L - 3]);
$dir_res = "$dir/resources/data-raw/Transcript-ids-cleaned/";
$dir_out = "$dir/resources/data-raw/protein-sequences-filtered/";

#opendir(DH, $dir_res);
#my @dirs = readdir(DH);
#close(DH);

# ------------------------
# Here we open the resource folder containing 
# all the species transcript ids
#
my $filecount = 0;
#my $folder_path = "$dir_res$folder";
opendir(DH, $dir_res);
my @spe_files = readdir(DH);
close(DH);

# ---------------------------
# Here we iterate over all the files in each 
# species folder
#
foreach my $file (@spe_files) {
    print "Currently downloading sequences for $file\n";
    # ---------------------------
    # this opens a file for each species to write the retrieved 
    # data to.
    #
    my $spe_output = "$dir_out$file.txt";
    open(FH, ">$spe_output");

    my $filepath = "$folder_path/$file";
    if (index($filepath, ".txt") != -1){
        $filecount += 1;
        print $filepath . "  ($filecount/23 files)\n";
        open(my $fh, '<:encoding(UTF-8)', $filepath) or die "Could not open file '$filepath' $!";

        # -------------------------------
        # then we open each file and iterate over every line in the 
        # file and use the second id (protein id) to query
        # the ensembl compara database.
        #
        while (my $row = <$fh>) {
            chomp $row;
            my @L = split ' ', $row;
            my $ID = @L[1];
            
            if (index($ID, "ENS") != -1){
                #fetch a Member
                my $seqmember = $seqmember_adaptor->fetch_by_stable_id($ID);
                print FH ">$ID\n";
                print FH $seqmember->sequence(),"*\n";                
            }
        }
    }   
}
close(FH);
