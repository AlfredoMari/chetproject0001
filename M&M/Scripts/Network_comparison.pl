use strict;
use Cwd;


# Alfredo Mari

##### GENERAL
# This script is meant to be used for comparing two Netwoks, i.e. two Edge tables which are the result of Matt's script Co_occurrence_Network_OTU_wPr.pl,
# recoverable in the same directory. The script picks up from the single Edge table the adj_rsq and from this builds up a distance matrix. All this is achieved through R.
# Afterwards the two distance matrices are compared through compare_distance_matrices.py which compares them and performs a Mantel correlogram.
# Mantel correlogram is giving you a score included between -1 and 1, -1 means strongly negative correlation, 1 means strong positive correlation, 0 no correlation.
# The script provides also a graph in which single distance classes are analyzed, and for each of them p value is calculated,
# this way is also possible to check wether the correlation is linear or not.
# Improvement compared to the previous versions:
# - Here is also calculated a general Mantel test to compare the overall matrices
# - I included the link to the backtracker script, able to track back the id/taxonomy of the interactors involved
# - I included a reverse test, basically adding to the N1 vs N2 comparison, also the N2 vs N1 comparison, this is not changing
# - I tried to unique all the interaction summary file for a better viewing and less confusion ###TO BE VERIFIED IF CORRECT###
#   the shape of the correlogram, but since the class indexes and the interactions displayed in the corr_results report belong onl to one
#   distance matrix, adding the reverse test we also add the other half of interactions which otherwise will be lost
# - I adapted it for Co_Net outputs
# - one can choose if to compare the networks through edge betweenness or weight, just select EDGEBET or WEIGHT as last passthrough for the commandline

##### USAGE:
# perl /projects/dep_psl/grp_kemen/programs/perl_scripts/Alfredo/Network_comparison_adj_rsq_debugged.pl 1STNETWORKNAME 1STNETWORK 2NDNETWORKNAME 2NDNETWORK
# Where:
# 1STNETWORKNAME: give a name to the first network you want to compare
# 1STNETWORK: /path/to/the/1stEdgetable.txt
# 2NDNETWORKNAME: give a name to the second network you are comparing with
# 2NDNETWORK: /path/to/the/2ndEdgetable.txt
# VARIABLE: it is either edgebetweenness or weight. values allowed are EDGEBET for edgebetweenness or WEIGHT for weight

my $fstNetName = shift;             
my $fstNet = shift;
my $sndNetName = shift;
my $sndNet = shift;
my $var = shift;
my $bsub = ". /opt/lsf/conf/profile.lsf";
my $Rc;
#Runs the extractor to generate the files out of the .csv table that cytoscape produces
my $out1 = "Extraction.out";
my $Ra = "bsub -o $out1 Rscript /biodata/dep_psl/grp_kemen/programs/perl_scripts/Alfredo/R_dependencies/Extractor.R $fstNet $sndNet";
my $Run = system("$bsub && $Ra") == 0 or die "Can't extract the networks";
sleep 1 until -e $out1;
my $fstNet1 = "minus_edge_table_div.txt";
my $sndNet1 = "plus_edge_table_div.txt";
#Calls R function in order to make a distance matrix based only on one value from the edgetable
my $out = "matrix_maker.out";
if ($var eq "WEIGHT") {
    $Rc = "bsub -o $out Rscript /biodata/dep_psl/grp_kemen/programs/perl_scripts/Alfredo/R_dependencies/matrix_maker_weight.R $fstNetName $fstNet1 $sndNetName $sndNet1";
}elsif($var eq "EDGEBET"){
    $Rc = "bsub -o $out Rscript /biodata/dep_psl/grp_kemen/programs/perl_scripts/Alfredo/R_dependencies/matrix_maker.R $fstNetName $fstNet1 $sndNetName $sndNet1";
}
    #execute it
my $Run = system("$bsub && $Rc") == 0 or die "Can't calculate the matrices";
sleep 1 until -e $out;
#Now we run the mantel correlogram to compare the two distance matrices 1 vs 2

my $fstMatrix = $fstNetName."_hollow_dist.txt";
my $sndMatrix = $sndNetName."_hollow_dist.txt";

my $fdir = "Comp_".$fstNetName."_vs_".$sndNetName;
my $sdir = "Comp_".$sndNetName."_vs_".$fstNetName;

my $qiime = ". /opt/share/software/packages/qiime-1.8.0/Activation_scripts/parallel/mpipzQIIME_1.8.0_new_parallel.s";
#1_vs_2
my $Mantel_corr_out1 = "Mantel_corr_".$fstNetName."_vs_".$sndNetName.".out";
my $mantel_corr1 = "bsub -o $Mantel_corr_out1 compare_distance_matrices.py --method mantel_corr -i $fstMatrix,$sndMatrix -o ".$fdir."/Mantel_correlogram/ -n 500";
my $Mantel_out1 = "Mantel_gen_".$fstNetName."_vs_".$sndNetName.".out";
my $Mantel1 = "bsub -o $Mantel_out1 compare_distance_matrices.py --method mantel -i $fstMatrix,$sndMatrix -o ".$fdir."/Mantel_general/ -n 500";

#2_vs_1
my $Mantel_corr_out2 = "Mantel_corr_".$sndNetName."_vs_".$fstNetName.".out";
my $mantel_corr2 = "bsub -o $Mantel_corr_out2 compare_distance_matrices.py --method mantel_corr -i $sndMatrix,$fstMatrix -o ".$sdir."/Mantel_correlogram/ -n 500";
my $Mantel_out2 = "Mantel_gen_".$sndNetName."_vs_".$fstNetName.".out";
my $Mantel2 = "bsub -o $Mantel_out2 compare_distance_matrices.py --method mantel -i $sndMatrix,$fstMatrix -o ".$sdir."/Mantel_general/ -n 500";

my $cmd0 = system('bash','-c',"$bsub && $qiime && $Mantel1 && $Mantel2") == 0
    or die "system failed: $?";


my $cmd1 = system('bash','-c',"$bsub && $qiime && $mantel_corr1 && $mantel_corr2") == 0
    or die "system failed: $?";
    
sleep 0 until -e $Mantel_out1;
sleep 0 until -e $Mantel_corr_out1;
sleep 0 until -e $Mantel_out2;
sleep 0 until -e $Mantel_corr_out2;

#Now we track back the interactions using the backtracker scripts, this will allow us
#to know which interactions do belong to which class, first 1 vs 2
chdir ("./".$fdir."/");
my $Backout1 = "Interaction_tracking_".$fstNetName."_vs_".$sndNetName.".out";
my $corresults1 = "./Mantel_correlogram/mantel_correlogram_results.txt";
my $Mat1 = "./../$fstMatrix";
my $Mat2 = "./../$sndMatrix";
my $backtracker1 = "bsub -o $Backout1 Rscript /biodata/dep_psl/grp_kemen/programs/perl_scripts/Alfredo/R_dependencies/Backtracker.R $fstNetName $Mat1 $sndNetName $Mat2 $corresults1";
my $cmd2 = system("$bsub && $backtracker1");
sleep 0 until -e $Backout1;

#Now let's do it for the reverse, 2 vs 1
chdir ("./../".$sdir."/");
my $Backout2 = "Interaction_tracking_".$sndNetName."_vs_".$fstNetName.".out";
my $corresults2 = "./Mantel_correlogram/mantel_correlogram_results.txt";
my $backtracker2 = "bsub -o $Backout2 Rscript /biodata/dep_psl/grp_kemen/programs/perl_scripts/Alfredo/R_dependencies/Backtracker.R  $sndNetName $Mat2 $fstNetName $Mat1 $corresults2";
my $cmd3 = system("$bsub && $backtracker2");


sleep 0 until -e $Backout2;

chdir ("./../");

# Now we open the interaction files conaining each one half of the interactions we are interested in
# and we merge them ito a unique file

my $Itrack1 = $fdir."/Interaction_description_".$fstNetName."_vs_".$sndNetName.".txt";
my $Itrack2 = $sdir."/Interaction_description_".$sndNetName."_vs_".$fstNetName.".txt";

my $out2= "merge_both.out";
my $Rcomm = "bsub -o $out2 Rscript /biodata/dep_psl/grp_kemen/programs/perl_scripts/Alfredo/R_dependencies/zweite.R $Itrack1 $Itrack2";

#execute it
my $run2 = system("$bsub && $Rcomm") ==0 or die "Can't merge the tw sets";
sleep 1 until -e $out2;
exit;
