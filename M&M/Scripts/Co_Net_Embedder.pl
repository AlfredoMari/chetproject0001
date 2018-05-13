#Alfredo Mari
#!/usr/bin/perl
use strict;
#This script is meant to embed all the required steps of CO_NET in order to run it on the command line.
#for more info on how does conet work, you can visit this page and follow the tutorial, or see the user friendly app on cytoscape.
#The reason why you may need to use this script is that, the GUI version of CONET is running only locally on your laptop and has restricted memory
#availability, meaning that a medium size network, to run it locally can take 6 to 8 hours. Plus, the GUI interface allows only one network at a time.
#this script allows you to use the cluster power to use coNet, meaning that the same network, to be calculated will take 20-30 min.
#It uses LSF, and this means you can submit multiple networks at a time, so in 1h you can have up to 10 Nets calculated.
#I queued the jobs in ioheavy, this make the whole procedure faster.
#Another improval is that the GUI versions asks you, for each network, to start three runs, one for the loading, one for the permutations, and one for the bootstrap.
#here all this passages are embedded sequentially, meaning that you just have to run the script once per network.

#prerun:
# export LIB=/biodata/dep_psl/grp_kemen/programs/cytoscape-unix-3.5.1/apps/CoNet3/lib
# export CLASSPATH=${CLASSPATH}:${LIB}/CoNet.jar
#V2 does not renormalize since re-normalization created problems
my $bsub = ". /opt/lsf/conf/profile.lsf";
my $input= shift;
my $features= shift;
my $guessparam= shift;   #suggested not to go below 1000, rare cases  500
my $filterparam= shift;   #suggested max 50, don't go below 10

my @split = split(/\./, $input);
my $name = $split[0];
my @splitname = split(/\_/, $name);
my $comp = $splitname[7];
#my $feat = $splitname[10]."_".$splitname[11];
my $tag = $comp."_gp_".$guessparam."_fp_".$filterparam;
my $dir = "./".$tag."/";
my $thresholdfile = $dir."thresholds.txt";
my $outnet= $dir."Network.gdl";
my $permnet= $dir."Network_permutated.gdl";
my $bootnet= $dir."Network_bootstraped.gdl";
my $permfile = $dir."Permutation_file.txt";
my $bootfile = $dir."Bootstrap_randomization_file.txt";
my $outdir = "./Outdir_".$tag."/";


my $diro = $tag."_dirmaker.out";
my $dircmd = "bsub -o $diro mkdir $outdir";
my $dirout = $outdir."dir.out";
my $dircmd1 = "bsub -o $dirout mkdir $dir";
#let's make directories
my $run = system('bash', '-c', "$bsub && $dircmd") == 0 or die "Can't create outdir"."\n";
sleep 1 until -e $diro;

my $run1 = system('bash','-c',"$bsub && $dircmd1") == 0 or die "Can't create directories"."\n";
sleep 1 until -e $dirout;

#infer the network
my $plout = $outdir.$tag."_Pre_load.out";
my $loadout = $outdir.$tag."_Load_net.out";
#new
my $preload = "bsub -q ioheavy -o $plout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --min 1 --inputformat biomtable --features $features --measure1 conf --multicorr none --minsupport 2 --measure2 supp --iterations 100 --correlnonrandp none --inference mrnet --topbottom --matrixtype count --stand col_norm --thresholdguessing edgeNumber --guessingparam $guessparam --resamplemethod shuffle_rows --ensemblemethods correl_pearson/correl_spearman/correl_kendall/dist_bray/dist_kullbackleibler --max 3 --format gdl --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --multigraph --higherleveltaxa --transposefeatures --matchfeaturesamples --keepfilteredrows --filter row_minocc/noinclusivetaxalinks --filterparameter $filterparam  --output $thresholdfile";
# old my $preload = "bsub -q ioheavy -o $plout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --matrixtype count --ensemblemethods correl_pearson/correl_spearman/correl_kendall/sim_mutInfo/dist_bray/dist_kullbackleibler --guessingparam $guessparam --thresholdguessing edgeNumber --resamplemethod shuffle_rows --max 3 --format gdl --stand col_norm --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --inputformat biomtable --features $features --correlnonrandp none --multicorr none --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --topbottom --multigraph --higherleveltaxa --transposefeatures --keepfilteredrows --matchfeaturesamples --filter row_minocc/noinclusivetaxalinks --filterparameter $filterparam --output $thresholdfile";
#my $load = "bsub -q ioheavy -o $loadout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --matrixtype count --randroutine none --ensemblemethods correl_pearson/correl_spearman/correl_kendall/sim_mutInfo/dist_bray/dist_kullbackleibler --resamplemethod shuffle_rows --max 3 --format gdl --stand col_norm --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --inputformat biomtable --features $features --correlnonrandp none --multicorr none --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --multigraph --higherleveltaxa --transposefeatures --keepfilteredrows --matchfeaturesamples --filter row_minocc/noinclusivetaxalinks --filterparameter $filterparam --output $outnet --ensembleparamfile $thresholdfile";
my $load = "bsub -q ioheavy -o $loadout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --min 1 --inputformat biomtable --features $features --measure1 conf --multicorr none --minsupport 2 --measure2 supp --iterations 100 --correlnonrandp none --inference mrnet --randroutine none --matrixtype count --stand col_norm --resamplemethod shuffle_rows --ensemblemethods correl_pearson/correl_spearman/correl_kendall/dist_bray/dist_kullbackleibler --max 3 --format gdl --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --multigraph --higherleveltaxa --transposefeatures --matchfeaturesamples --keepfilteredrows --filter row_minocc/noinclusivetaxalinks --filterparameter $filterparam --output $outnet --ensembleparamfile $thresholdfile";

my $run2 = system('bash','-c',"$bsub && $preload") == 0 or die "Can't infer the network, preload, line 62"."\n";
sleep 1 until -e $plout;
my $run3 = system ('bash','-c',"$bsub && $load") == 0 or die "Can't infer the network, load, line 63"."\n";
sleep 1 until -e $loadout;

#pemutate the network
my $perm_plout = $outdir.$tag."_Perm_Preload.out";
my $perm_loadout = $outdir.$tag."_Perm_Load_net.out";
#old my $perm_preload = "bsub -q ioheavy -o $perm_plout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --matrixtype count --ensemblemethods correl_pearson/correl_spearman/correl_kendall/sim_mutInfo/dist_bray/dist_kullbackleibler --guessingparam $guessparam --thresholdguessing edgeNumber --resamplemethod shuffle_rows --max 3 --format gdl --stand col_norm --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --inputformat biomtable --features $features --correlnonrandp none --multicorr none --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --topbottom --multigraph --renorm --higherleveltaxa --transposefeatures --keepfilteredrows --matchfeaturesamples --filter row_minocc/noinclusivetaxalinks/rand --filterparameter $filterparam  --output $thresholdfile";
# renormalized my $perm_preload = "bsub -q ioheavy -o $perm_plout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --min 1 --inputformat biomtable --features $features --measure1 conf --multicorr none --minsupport 2 --measure2 supp --iterations 100 --correlnonrandp none --inference mrnet --topbottom --matrixtype count --stand col_norm --thresholdguessing edgeNumber --guessingparam $guessparam --resamplemethod shuffle_rows --ensemblemethods correl_pearson/correl_spearman/correl_kendall/dist_bray/dist_kullbackleibler --max 3 --format gdl --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --multigraph --higherleveltaxa --renorm --transposefeatures --matchfeaturesamples --keepfilteredrows --filter row_minocc/noinclusivetaxalinks/rand --filterparameter $filterparam  --output $thresholdfile";
# old my $perm_load = "bsub -q ioheavy -o $perm_loadout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --randscorefile $permfile --matrixtype count --randroutine edgeScores --ensemblemethods correl_pearson/correl_spearman/correl_kendall/sim_mutInfo/dist_bray/dist_kullbackleibler --resamplemethod shuffle_rows --max 3 --format gdl --stand col_norm --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --inputformat biomtable --features $features --correlnonrandp none --multicorr none --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --multigraph --renorm --scoreexport --higherleveltaxa --transposefeatures --keepfilteredrows --matchfeaturesamples --filter row_minocc/noinclusivetaxalinks/rand --filterparameter $filterparam --output $permnet --ensembleparamfile $thresholdfile";
my $perm_preload = "bsub -q ioheavy -o $perm_plout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --thresholdguessing edgeNumber --guessingparam $guessparam --resamplemethod shuffle_rows --scoremergestrategy mean --topbottom --format gdl --ensemblemethods correl_pearson/correl_spearman/dist_bray/dist_kullbackleibler --max 3 --stand col_norm --matrixtype count --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --inputformat biomtable --features $features --correlnonrandp none --multicorr none --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --keepfilteredrows --matchfeaturesamples --transposefeatures --higherleveltaxa --multigraph --filter row_minocc/noinclusivetaxalinks/rand --filterparameter $filterparam  --output $thresholdfile";
# renormalized my $perm_load = "bsub -q ioheavy -o $perm_loadout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --min 1 --inputformat biomtable --features $features --measure1 conf --multicorr none --minsupport 2 --measure2 supp --iterations 100 --correlnonrandp none --inference mrnet --randroutine edgeScores --randscorefile $permfile --matrixtype count --stand col_norm --resamplemethod shuffle_rows --ensemblemethods correl_pearson/correl_spearman/correl_kendall/dist_bray/dist_kullbackleibler --max 3 --format gdl --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --multigraph --higherleveltaxa --renorm --scoreexport --transposefeatures --matchfeaturesamples --keepfilteredrows --filter row_minocc/noinclusivetaxalinks/rand --filterparameter $filterparam --output $permnet --ensembleparamfile $thresholdfile";
my $perm_load = "bsub -q ioheavy -o $perm_loadout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --resamplemethod shuffle_rows --scoremergestrategy mean --format gdl --ensemblemethods correl_pearson/correl_spearman/dist_bray/dist_kullbackleibler --max 3 --randscorefile $permfile --randroutine edgeScores --stand col_norm --matrixtype count --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --inputformat biomtable --features $features --correlnonrandp none --multicorr none --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --keepfilteredrows --matchfeaturesamples --transposefeatures --scoreexport --higherleveltaxa --multigraph --filter row_minocc/noinclusivetaxalinks/rand --filterparameter $filterparam --output $permnet --ensembleparamfile $thresholdfile";
my $run4 = system('bash','-c',"$bsub && $perm_preload") == 0 or die "Can't permutate the network, preload, line 77"."\n";
sleep 1 until -e $perm_plout;

my $run5 = system('bash','-c',"$bsub && $perm_load") == 0 or die "Can't permutate the network, load, line 78"."\n";
sleep 1 until -e $perm_loadout;

#bootstrap the network
my $boot_plout = $outdir.$tag."_Boot_Preload.out";
my $boot_loadout = $outdir.$tag."_Boot_Load_net.out";
#old my $boot_preload = "bsub -q ioheavy -o $boot_plout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --matrixtype count --ensemblemethods correl_pearson/correl_spearman/correl_kendall/sim_mutInfo/dist_bray/dist_kullbackleibler --nulldistribfile $permfile --guessingparam $guessparam --thresholdguessing edgeNumber --resamplemethod bootstrap --max 3 --format gdl --stand col_norm --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --inputformat biomtable --features $features --correlnonrandp none --multicorr benjaminihochberg --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --pvaluemerge brown --topbottom --multigraph --higherleveltaxa --transposefeatures --keepfilteredrows --matchfeaturesamples --filter row_minocc/noinclusivetaxalinks/rand/confidence_boot --filterparameter $filterparam  --output $thresholdfile";
my $boot_preload = "bsub -q ioheavy -o $boot_plout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --min 1 --inputformat biomtable --features $features --measure1 conf --multicorr benjaminihochberg --minsupport 2 --measure2 supp --iterations 100 --correlnonrandp none --inference mrnet --pvaluemerge brown --topbottom --matrixtype count --stand col_norm --thresholdguessing edgeNumber --guessingparam $guessparam --resamplemethod bootstrap --ensemblemethods correl_pearson/correl_spearman/correl_kendall/dist_bray/dist_kullbackleibler --max 3 --format gdl --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --multigraph --higherleveltaxa --transposefeatures --matchfeaturesamples --keepfilteredrows --filter row_minocc/noinclusivetaxalinks/rand/confidence_boot --filterparameter $filterparam  --output $thresholdfile";
#old my $boot_load = "bsub -q ioheavy -o $boot_plout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --randscorefile $bootfile --matrixtype count --randroutine edgeScores --ensemblemethods correl_pearson/correl_spearman/correl_kendall/sim_mutInfo/dist_bray/dist_kullbackleibler --nulldistribfile $permfile --resamplemethod bootstrap --max 3 --format gdl --stand col_norm --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --inputformat biomtable --features $features --correlnonrandp none --multicorr benjaminihochberg --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --pvaluemerge brown --multigraph --scoreexport --higherleveltaxa --transposefeatures --keepfilteredrows --matchfeaturesamples --filter row_minocc/noinclusivetaxalinks/rand/confidence_boot --filterparameter $filterparam --output $bootnet --ensembleparamfile $thresholdfile";
my $boot_load = "bsub -q ioheavy -o $boot_loadout java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $input --scoremergestrategy mean --min 1 --inputformat biomtable --features $features --measure1 conf --multicorr benjaminihochberg --minsupport 2 --measure2 supp --iterations 100 --correlnonrandp none --inference mrnet --pvaluemerge brown --randroutine edgeScores --randscorefile $bootfile --matrixtype count --stand col_norm --resamplemethod bootstrap --ensemblemethods correl_pearson/correl_spearman/correl_kendall/dist_bray/dist_kullbackleibler --max 3 --format gdl --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --multigraph --higherleveltaxa --scoreexport --transposefeatures --matchfeaturesamples --keepfilteredrows --filter row_minocc/noinclusivetaxalinks/rand/confidence_boot --filterparameter $filterparam --output $bootnet --ensembleparamfile $thresholdfile";

my $run6 = system('bash','-c',"$bsub && $boot_preload") == 0 or die "Can't bootstrap the network, preload, line 89"."\n";
sleep 1 until -e $boot_plout;

my $run7 = system('bash','-c',"$bsub && $boot_load") == 0 or die "Can't bootstrap the network, load, line 90"."\n";
sleep 1 until -e $boot_loadout;



exit;



