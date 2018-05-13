
use strict;

#should be run from the directory from which one picks the otus (usearch_BacV3, for example)

my $primer = shift;
my $otumap = "./pick_otus_all_prefiltered/final_otu_map_fixed_mc2.txt";
my $tax = "./RDP_assigned_taxa/final_rep_set_fixed_tax_assignments.txt";
my $otutab = $primer."_OTU_table_general_native_pr_mc2.biom";
my $qiime = ". /opt/share/software/packages/qiime-1.8.0/Activation_scripts/parallel/mpipzQIIME_1.8.0_new_parallel.s ";
my $bsub = ". /opt/lsf/conf/profile.lsf ";

#first we create an outdir where to store the outfiles
my $mainoutdir = "./Outdir/";
my $mainout = "outdir.out";
my $cmd = "bsub -o $mainout mkdir $mainoutdir";
my $rune = system ('bash','-c',"$qiime && $bsub && $cmd") == 0 or die "Can't create Outdir";
sleep 1 until -e $mainout;

my $tabout = $mainoutdir."make_otu_tab_".$primer.".out";
my $otumaker = "bsub -o $tabout make_otu_table.py -i $otumap -t $tax -o $otutab";

my $run = system ('bash','-c',"$qiime && $bsub && $otumaker") == 0 or die "Can't make the otu table"."\n";

sleep 1 until -e $tabout;

#now we filter out the otu table based on the occurrence of the otus in the samples

my $filtab= $primer."_OTU_table_general_native_pr_mc2_s2_n50.biom";
my $filtout = $mainoutdir.$primer."general_filter.out";
my $gen_filter = "bsub -o $filtout filter_otus_from_otu_table.py -i $otutab -o $filtab -s 2 -n 50";

my $run1 = system ('bash','-c',"$bsub && $qiime && $gen_filter") == 0 or die "Can't perform the general filter [s2 n50]"."\n";

sleep 1 until -e $filtout;

#now we keep only the otus belonging to the samples that are present in the mapfile

my $filtmap = $primer."_OTU_table_general_core_pr_mc2_s2_n50.biom";
my $filtmapout = $mainoutdir.$primer."_mapfile_core_filter.out";
my $map_filter = "bsub -o $filtmapout filter_samples_from_otu_table.py -i $filtab -o $filtmap --sample_id_fp /netscratch/dep_psl/grp_kemen/netscratch_alf/Master_map_file__NL_PCOA.txt";

my $run2 = system ('bash','-c',"$bsub && $qiime && $map_filter") == 0 or die "Can't perform the mapfile trimming"."\n";

sleep 1 until -e $filtmapout;

#now we create the endo and wl table (and also the root endo and epi, we will merge them later on anyway)

my $Endotab = $primer."_OTU_table_general_core_pr_mc2_s2_n50_Endo.biom";
my $Wltab = $primer."_OTU_table_general_core_pr_mc2_s2_n50_Wl.biom";
my $REndotab =$primer."_OTU_table_general_core_pr_mc2_s2_n50_REndo.biom";
my $REpitab =$primer."_OTU_table_general_core_pr_mc2_s2_n50_REpi.biom";

my $Endout = $mainoutdir.$primer."_Endotab.out";
my $Wlout = $mainoutdir.$primer."_Wltab.out";
my $REndout = $mainoutdir.$primer."_REndotab.out";
my $REpiout = $mainoutdir.$primer."_REpitab.out";

my $Endomaker= "bsub -o $Endout filter_samples_from_otu_table.py -i $filtmap -o $Endotab -m /netscratch/dep_psl/grp_kemen/netscratch_alf/Master_map_file__NL_PCOA.txt -s Compartment:Endo";
my $Wlmaker="bsub -o $Wlout filter_samples_from_otu_table.py -i $filtmap -o $Wltab -m /netscratch/dep_psl/grp_kemen/netscratch_alf/Master_map_file__NL_PCOA.txt -s Compartment:Wl";
my $REndomaker="bsub -o $REndout filter_samples_from_otu_table.py -i $filtmap -o $REndotab -m /netscratch/dep_psl/grp_kemen/netscratch_alf/Master_map_file__NL_PCOA.txt -s Compartment:R_Endo";
my $REpimaker="bsub -o $REpiout filter_samples_from_otu_table.py -i $filtmap -o $REpitab -m /netscratch/dep_psl/grp_kemen/netscratch_alf/Master_map_file__NL_PCOA.txt -s Compartment:R_Epi";

my $run3 = system ('bash','-c',"$bsub && $qiime && $Endomaker && $Wlmaker && $REndomaker && $REpimaker") == 0 or die "Can't divide the tables by compartment tables"."\n";

sleep 1 until -e $Endout;
sleep 1 until -e $Wlout;
sleep 1 until -e $REndout;
sleep 1 until -e $REpiout;

#now we run the otu pipe script in order to have the epiphytic table of the leaf ad also to have the summaries and the tax filter
my $OTUout = $mainoutdir."OTU_pipe.out";
my $OTUpipe = "bsub -o $OTUout perl /projects/dep_psl/grp_kemen/programs/perl_scripts/Alfredo/OTU_filter_pipe_mod.pl $Endotab $Wltab $REndotab $REpitab";
my $run4 = system ('bash','-c',"$bsub && $OTUpipe") == 0 or die "Can't perform OTU pipe.pl "."\n";

sleep 1 until -e $OTUout;


exit;


