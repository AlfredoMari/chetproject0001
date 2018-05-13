#Alfredo Mari
#This script is made for Filter the provided OTU table in order to prune the Epiphytic organisms in case one has only Whole leaf and Endophytes otu tables.
#The embedded r scripts is processing the OTU table not summarised, it produces three OTU tables, one Endo, one Epi, one Wl all wothout 0s and artifacts otus.
#Going forward, the script cleans up from the unwanted taxa, like arabidopsis, mitochondria DNA, or Fungi otus for 18S.
#The scripts proceeds in finding automatically a rarefaction threshold via before performaing a summary of the reads per sample, (accessible in the directory Summaries/)
#the threshold is chosen via scanning the first seven samples, if none of them is above 500 reads per sample (300 for epiphytes), 500 reads (300 for epiphytes) are chosen by default,
#if instead one of the seven is above 500 (300 for epiphytes) reads, that number is chosen as a threshold.
#This way we are sure not trash more than seven samples, unless the default choice is used.

###IMPORTANT####
#The summaries directory anyhow should be checked, since each data set is different and seven samples might be not enough as a chek basis.
#afterwards, the script is summarising the otus by taxa, and putting them in the Taxa_summarised directory.

###NOTE###
#that the script is taking in account the different format of the taxa assignement for protists and for others primer set, so can be used also for them,
#just keep in mind that level 4 for protists means Class and not Family.
#as last, the script is merging the epi and Endo OTU table (the rarefied, not summarised) in only one otu table, that can be used for R plots and PCOA constraint.

###### REQUIREMENTS:
#please run the following commands for running the script:
# export PERL5LIB=$PERL5LIB:/projects/dep_psl/grp_kemen/programs/PerlMods/Statistics-R-0.34/lib/:/projects/dep_psl/grp_kemen/programs/PerlMods/Regexp-Common-2013031301/lib/:/projects/dep_psl/grp_kemen/programs/PerlMods/
# export PATH=/projects/dep_psl/grp_kemen/programs/R-3.3.1/bin:$PATH
# export R_LIBS="/projects/dep_psl/grp_kemen/programs/myRLibs/"

#in order to run everything one should have ready:
#a OTU table containing Whole leaf samples, WLTABLE
#a OTU table containing Endophytic camples, ENDOTABLE
#both of them have to be in .biom format

############# !!!!!!!!!! VERY IMPORTANT !!!!!!!!!! #############

# Do not source qiime before running the script, otherwise the whole thing won't work!!

############# !!!!!!!!!! VERY IMPORTANT !!!!!!!!!! #############

#USAGE

# perl /projects/dep_psl/grp_kemen/programs/perl_scripts/Alfredo/OTU_filter_pipe.pl ENDOTABLE WLTABLE

#script starts:
# set packages 
use strict;
use Cwd;
use Statistics::R;
use List::Util qw (min);
my $R = Statistics::R->new();
my $bsub = ". /opt/lsf/conf/profile.lsf";
my $qiime= ". /opt/share/software/packages/qiime-1.8.0/Activation_scripts/parallel/mpipzQIIME_1.8.0_new_parallel.s";
#custom input
my $Endotable=shift;
my $Wltable=shift;
my $REndotab = shift;
my $REpitab = shift;


#define names of processing files and outputs
my @split = split(/\_/, $Endotable);
my $primer = $split[0];
my @splitEndo = split(/\./, $Endotable);
my @splitWl = split(/\./, $Wltable);

my $Endoname = $splitEndo[0];
my $Wlname = $splitWl[0];

#define output directories
my $Maindir = "./$primer/";
my $NoRardir = "./$primer/Not_rarefied/";
#my $Rardir = "./$primer/Rarefied/";
my $Sumdir = "./$primer/Summaries/";
#my $Taxasumdir = "./$primer/Rarefied/Taxa_summarised/";
my $outdir = "./".$primer."_outdir/";


my $Endotxt = $NoRardir.$Endoname.".txt";
my $Wltxt = $NoRardir.$Wlname.".txt";

my $FinalEndo_txt = $NoRardir.$primer."_native_NO_0_Endo.txt";
my $FinalEpi_txt = $NoRardir.$primer."_native_NO_0_Epi.txt";
my $FinalWl_txt = $NoRardir.$primer."_native_NO_0_Wl.txt";

my $Endo_partial_biom = $NoRardir.$primer."_native_NO_0_Endo.biom";
my $Epi_partial_biom = $NoRardir.$primer."_native_NO_0_Epi.biom";

my $FinalEndo_biom = $NoRardir.$primer."_native_Rincluded_NO0_Endo.biom";
my $FinalEpi_biom = $NoRardir.$primer."_native_Rincluded_NO0_Epi.biom";
my $FinalWl_biom = $NoRardir.$primer."_native_NO_0_Wl.biom";

my $EndoSummary = $Sumdir.$primer."_native_NO_0_Endo_taxacleaned_summary.txt";
my $WlSummary = $Sumdir.$primer."_native_NO_0_Wl_taxacleaned_summary.txt";
my $EpiSummary = $Sumdir.$primer."_native_NO_0_Epi_taxacleaned_summary.txt";

my $TaxacleanEndo_b = $NoRardir.$primer."_NO_0_taxa_cleaned_Endo.biom";
my $TaxacleanWl_b = $NoRardir.$primer."_NO_0_taxa_cleaned_Wl.biom";
my $TaxacleanEpi_b = $NoRardir.$primer."_NO_0_taxa_cleaned_Epi.biom";


#my $TaxacleanEndo_t = $NoRardir.$primer."_NO_0_taxa_cleaned_Endo.txt";
#my $TaxacleanWl_t = $NoRardir.$primer."_NO_0_taxa_cleaned_Wl.txt";
#my $TaxacleanEpi_t = $NoRardir.$primer."_NO_0_taxa_cleaned_Epi.txt";

#let's start with the outdir

my $dashout = "Outdir_".$primer.".out";
my $out = "bsub -o $dashout mkdir $outdir";

my $run0 = system('bash','-c', "$bsub && $out") == 0
  or die "Can't create outdir, see line 100"."\n";
sleep 1 until -e $dashout;

my $dirout = $outdir."Maindir.out";
my $dirout1 = $outdir."NoRardir.out";
#my $dirout2 = $outdir."Rardir.out";
#my $dirout3 = $outdir."Sumdir.out";
#my $dirout4 = $outdir."Taxasumdir.out";

my $dirmaker =  "bsub -o $dirout mkdir $Maindir";
my $dirmaker1 = "bsub -o $dirout1 mkdir $NoRardir";
#my $dirmaker2 = "bsub -o $dirout2 mkdir $Rardir";
#my $dirmaker3 = "bsub -o $dirout3 mkdir $Sumdir";
#my $dirmaker4 = "bsub -o $dirout4 mkdir $Taxasumdir";

my $run1 = system('bash', '-c', "$bsub && $dirmaker") == 0
  or die "Can't create Main directory, see line 116"."\n";

sleep 1 until -e $dirout;
my $run2 = system('bash', '-c', "$bsub && $dirmaker1") == 0
  or die "Can't create other directories, see line 120"."\n";

#sleep 1 until -e $dirout2;
#my $run2bis = system('bash', '-c', "$bsub && $dirmaker4") == 0
  or die "Can't create other directories, see line 124"."\n";

sleep 1 until -e $dirout1;
#sleep 1 until -e $dirout3;
#sleep 1 until -e $dirout4;

#first of, we need to convert the .biom tables in .txt tables:
my $tout1 = $outdir."1st_Txt_conv_raw_Endo.out";
my $tout2 = $outdir."1st_Txt_conv_raw_Wl.out";

my $txtconv1= "bsub -o $tout1 biom convert -b -i $Endotable -o $Endotxt --table-type \"otu table\" --header-key taxonomy";
my $txtconv2= "bsub -o $tout2 biom convert -b -i $Wltable -o $Wltxt --table-type \"otu table\" --header-key taxonomy";

my $run23 = system('bash', '-c', "$bsub && $qiime && $txtconv1 && $txtconv2") == 0
  or die "Can't convert biom to txt, see line 138"."\n";

sleep 1 until -e $tout1;
sleep 1 until -e $tout2;


my $Rcmd = "
#read in the files
dfwh <- read.csv(file = \"$Wltxt\", header =T, row.names =1, sep = \"\t\", dec = \".\", skip =1)
dfen <- read.csv(file = \"$Endotxt\", header =T, row.names =1, sep = \"\t\", dec = \".\", skip =1)

#consider only the dataframe without taxonomy
dfwh_t <- dfwh[,-which(names(dfwh)==\"taxonomy\")]
dfen_t <- dfen[,-which(names(dfen)==\"taxonomy\")]

ind_wh <- rowSums(dfwh_t) 
ind_en <- rowSums(dfen_t)

rev_wh <- which(ind_wh == 0, arr.ind = TRUE)
rev_en <- which(ind_en == 0, arr.ind = TRUE)

rev_wh <- as.vector(rev_wh)
rev_en <- as.vector(rev_en)

Wh_no_0 <- dfwh[-rev_wh,]
En_no_0 <- dfen[-rev_en,]

#now the trimming
ind_ep <- rownames(Wh_no_0) %in% rownames(En_no_0)
finder_ep <- which(ind_ep == \"FALSE\", arr.ind = TRUE)
Ep_no_0 <- Wh_no_0[finder_ep,]

#now formatting back to OTU table stile, so that we can transform all three OTU table back to biom
mat_en <- as.matrix(En_no_0)
Endophytes_end <- data.frame(\"OTU ID\" = rownames(mat_en),mat_en)
write.table(Endophytes_end, file = \"$FinalEndo_txt\", sep= \"\t\", quote =F,row.names = F, dec =\".\")

mat_wh <- as.matrix(Wh_no_0)
Whole_leaf_end <- data.frame(\"OTU ID\" = rownames(mat_wh),mat_wh)
write.table(Whole_leaf_end, file = \"$FinalWl_txt\", sep= \"\t\", quote =F,row.names = F, dec =\".\")

mat_ep <- as.matrix(Ep_no_0)
Epiphytes_end <- data.frame(\"OTU ID\" = rownames(mat_ep),mat_ep)
write.table(Epiphytes_end, file = \"$FinalEpi_txt\", sep= \"\t\", quote =F,row.names = F, dec =\".\")
";

$R -> send("$Rcmd");

sleep 1 until -e $FinalEndo_txt; 
sleep 1 until -e $FinalWl_txt;
sleep 1 until -e $FinalEpi_txt;
sleep 5;

#now let's create the biom file for all of them
my $b_out1 = $outdir."biom_conv_Native_Endo.out";
my $b_out2 = $outdir."biom_conv_Native_Wl.out";
my $b_out3 = $outdir."biom_conv_Native_Epi.out";

my $biomconv1 = "bsub -o $b_out1 biom convert -i $FinalEndo_txt -o $Endo_partial_biom --table-type \"OTU table\" --process-obs-metadata taxonomy";
my $biomconv2 = "bsub -o $b_out2 biom convert -i $FinalWl_txt -o $FinalWl_biom --table-type \"OTU table\" --process-obs-metadata taxonomy";
my $biomconv3 = "bsub -o $b_out3 biom convert -i $FinalEpi_txt -o $Epi_partial_biom --table-type \"OTU table\" --process-obs-metadata taxonomy";

my $run3 = system('bash', '-c', "$bsub && $qiime && $biomconv1 && $biomconv2 && $biomconv3") == 0
  or die "Can't convert the output table (Endo, Wl, Epi) back to .biom, see line 201"."\n";
sleep 1 until -e $b_out1; 
sleep 1 until -e $b_out2;
sleep 1 until -e $b_out3;
#Now we merge the root Endo with the leaf endo, and the epiphytes of the root with the new calculated epiphytes of the leaf.

my $mergendo_out = $outdir."Endophytes_merging.out";
my $mergepi_out = $outdir."Epiphytes_merging.out";

my $mergendo= "bsub -o $mergendo_out merge_otu_tables.py -i $Endo_partial_biom,$REndotab -o $FinalEndo_biom ";
my $mergepi= "bsub -o $mergepi_out merge_otu_tables.py -i $Epi_partial_biom,$REpitab -o $FinalEpi_biom ";

my $run_w = system ('bash','-c',"$bsub  && $qiime && $mergendo && $mergepi") == 0 or die "Can't merge the tables line 220";
sleep 1 until -e


#now let's work with the new formed OTU tables, they are native, so they need to be cleaned from the unwanted taxa
my $cleancmd1;
my $cleancmd2;
my $cleancmd3;


my $cleanout1 = $outdir."Filter_taxa_Endo.out";
my $cleanout2 = $outdir."Filter_taxa_Wl.out";
my $cleanout3 = $outdir."Filter_taxa_Epi.out";




if ($primer eq "BV3" || $primer eq "BV5"){

   $cleancmd1 ="bsub -o $cleanout1 filter_taxa_from_otu_table.py -i $FinalEndo_biom -o $TaxacleanEndo_b -n Unassigned,f__mitochondria,c__Chloroplast,o__Rickettsiales ";
   $cleancmd2 ="bsub -o $cleanout2 filter_taxa_from_otu_table.py -i $FinalWl_biom -o $TaxacleanWl_b -n Unassigned,f__mitochondria,c__Chloroplast,o__Rickettsiales ";
   $cleancmd3 ="bsub -o $cleanout3 filter_taxa_from_otu_table.py -i $FinalEpi_biom -o $TaxacleanEpi_b -n Unassigned,f__mitochondria,c__Chloroplast,o__Rickettsiales ";

}elsif ($primer eq "FITS2" || $primer eq "Ftrad"){

   $cleancmd1 ="bsub -o $cleanout1 filter_taxa_from_otu_table.py -i $FinalEndo_biom -o $TaxacleanEndo_b -n Unassigned,c__Oomycetes,k__Viridiplantae ";
   $cleancmd2 ="bsub -o $cleanout2 filter_taxa_from_otu_table.py -i $FinalWl_biom -o $TaxacleanWl_b -n Unassigned,c__Oomycetes,k__Viridiplantae ";
   $cleancmd3 ="bsub -o $cleanout3 filter_taxa_from_otu_table.py -i $FinalEpi_biom -o $TaxacleanEpi_b -n Unassigned,c__Oomycetes,k__Viridiplantae ";


}elsif ($primer eq "OITS2" || $primer eq "Otrad"){

   $cleancmd1 ="bsub -o $cleanout1 filter_taxa_from_otu_table.py -i $FinalEndo_biom -o $TaxacleanEndo_b -n Unassigned,k__Fungi ";
   $cleancmd2 ="bsub -o $cleanout2 filter_taxa_from_otu_table.py -i $FinalWl_biom -o $TaxacleanWl_b -n Unassigned,k__Fungi ";
   $cleancmd3 ="bsub -o $cleanout3 filter_taxa_from_otu_table.py -i $FinalEpi_biom -o $TaxacleanEpi_b -n Unassigned,k__Fungi ";

}elsif ($primer eq "PV4" || $primer eq "PV9"){
   $cleancmd1 ="bsub -o $cleanout1 filter_taxa_from_otu_table.py -i $FinalEndo_biom -o $TaxacleanEndo_b -n Unassigned,k__Fungi,k__Streptophyta,k__Stramenopiles_X ";
   $cleancmd2 ="bsub -o $cleanout2 filter_taxa_from_otu_table.py -i $FinalWl_biom -o $TaxacleanWl_b -n Unassigned,k__Fungi,k__Streptophyta,k__Stramenopiles_X ";
   $cleancmd3 ="bsub -o $cleanout3 filter_taxa_from_otu_table.py -i $FinalEpi_biom -o $TaxacleanEpi_b -n Unassigned,k__Fungi,k__Streptophyta,k__Stramenopiles_X ";

}
#Let's run it

my $run4 = system('bash', '-c', "$bsub && $qiime && $cleancmd1 && $cleancmd2 && $cleancmd3") == 0
  or die "Can't clean the taxa (Endo, Wl, Epi) see line 242"."\n";

sleep 1 until -e $cleanout1;
sleep 1 until -e $cleanout2;
sleep 1 until -e $cleanout3;


#Now we summarize, in order to understand how many reads per sample we have, and choose automatically the rarefactionn threshold
my $sumout1= $outdir."Summary_Endo.out";
my $sumout2= $outdir."Summary_Wl.out";
my $sumout3= $outdir."Summary_Epi.out";


my $getsum1 = "bsub -o $sumout1 biom summarize-table -i $TaxacleanEndo_b -o $EndoSummary";
my $getsum2 = "bsub -o $sumout2 biom summarize-table -i $TaxacleanWl_b -o $WlSummary";
my $getsum3 = "bsub -o $sumout3 biom summarize-table -i $TaxacleanEpi_b -o $EpiSummary";

my $run5 = system('bash', '-c', "$bsub && $qiime && $getsum1 && $getsum2 && $getsum3") == 0
  or die "Can't get the summary (Endo, Wl, Epi) see line 259"."\n";

sleep 1 until -e $sumout1;
sleep 1 until -e $sumout2;
sleep 1 until -e $sumout3;



##now we open the files and choose the threshold
#my $rardepth_w;
#open (SUMWL, $WlSummary) or die "Can't open the summary file Wl, see line 267"."\n";
#my @data_w = <SUMWL>;
#chomp @data_w;
#my $firstline_w = @data_w[16];
#my $secondline_w = @data_w[17];
#my $thirdline_w = @data_w[18];
#my $fourthline_w = @data_w[19];
#my $fifthline_w = @data_w[20];
#my $sixthline_w = @data_w[21];
#my $seventhline_w = @data_w[22];
#my @splitfirstline_w = split(/\:/,$firstline_w);
#my @splitsecondline_w = split(/\:/,$secondline_w);
#my @splitthirdline_w = split(/\:/,$thirdline_w);
#my @splitfourthline_w = split(/\:/,$fourthline_w);
#my @splitfifthline_w = split(/\:/,$fifthline_w);
#my @splitsixthline_w = split(/\:/,$sixthline_w);
#my @splitseventhline_w = split(/\:/,$seventhline_w);
#
#my $firstvalue_w = $splitfirstline_w[1];
#my $secondvalue_w = $splitsecondline_w[1];
#my $thirdvalue_w = $splitthirdline_w[1];
#my $fourthvalue_w = $splitfourthline_w[1];
#my $fifthvalue_w = $splitfifthline_w[1];
#my $sixthvalue_w = $splitsixthline_w[1];
#my $seventhvalue_w = $splitseventhline_w[1];
#
#my @splitvalue_w1 = split (/\./, $firstvalue_w);
#my @splitvalue_w2 = split (/\./, $secondvalue_w);
#my @splitvalue_w3 = split (/\./, $thirdvalue_w);
#my @splitvalue_w4 = split (/\./, $fourthvalue_w);
#my @splitvalue_w5 = split (/\./, $fifthvalue_w);
#my @splitvalue_w6 = split (/\./, $sixthvalue_w);
#my @splitvalue_w7 = split (/\./, $seventhvalue_w);
#
#my $firstdepth_w = $splitvalue_w1[0];
#my $seconddepth_w = $splitvalue_w2[0];
#my $thirddepth_w = $splitvalue_w3[0];
#my $fourthdepth_w = $splitvalue_w4[0];
#my $fifthdepth_w = $splitvalue_w5[0];
#my $sixthdepth_w = $splitvalue_w6[0];
#my $seventhdepth_w = $splitvalue_w7[0];
#
#$firstdepth_w =~ s/\s//;
#$seconddepth_w =~ s/\s//;
#$thirddepth_w =~ s/\s//;
#$fourthdepth_w =~ s/\s//;
#$fifthdepth_w =~ s/\s//;
#$sixthdepth_w =~ s/\s//;
#$seventhdepth_w =~ s/\s//;
#
#my $min_w = min grep {$_ > 500} $firstdepth_w, $seconddepth_w, $thirddepth_w, $fourthdepth_w, $fifthdepth_w, $sixthdepth_w, $seventhdepth_w;
#
#if (defined($min_w) && $min_w >= 500){
#  $rardepth_w = $min_w;
#}else{
#  $rardepth_w = 500;
#}
#close SUMWL;
#
#my $rardepth_ep;
#open (SUMEPI, $EpiSummary) or die "Can't open the summary file Epi, see line 327"."\n";
#my @data_ep = <SUMEPI>;
#chomp @data_ep;
#my $firstline_ep = @data_ep[16];
#my $secondline_ep = @data_ep[17];
#my $thirdline_ep = @data_ep[18];
#my $fourthline_ep = @data_ep[19];
#my $fifthline_ep = @data_ep[20];
#my $sixthline_ep = @data_ep[21];
#my $seventhline_ep = @data_ep[22];
#my @splitfirstline_ep = split(/\:/,$firstline_ep);
#my @splitsecondline_ep = split(/\:/,$secondline_ep);
#my @splitthirdline_ep = split(/\:/,$thirdline_ep);
#my @splitfourthline_ep = split(/\:/,$fourthline_ep);
#my @splitfifthline_ep = split(/\:/,$fifthline_ep);
#my @splitsixthline_ep = split(/\:/,$sixthline_ep);
#my @splitseventhline_ep = split(/\:/,$seventhline_ep);
#
#my $firstvalue_ep = $splitfirstline_ep[1];
#my $secondvalue_ep = $splitsecondline_ep[1];
#my $thirdvalue_ep = $splitthirdline_ep[1];
#my $fourthvalue_ep = $splitfourthline_ep[1];
#my $fifthvalue_ep = $splitfifthline_ep[1];
#my $sixthvalue_ep = $splitsixthline_ep[1];
#my $seventhvalue_ep = $splitseventhline_ep[1];
#
#my @splitvalue_ep1 = split (/\./, $firstvalue_ep);
#my @splitvalue_ep2 = split (/\./, $secondvalue_ep);
#my @splitvalue_ep3 = split (/\./, $thirdvalue_ep);
#my @splitvalue_ep4 = split (/\./, $fourthvalue_ep);
#my @splitvalue_ep5 = split (/\./, $fifthvalue_ep);
#my @splitvalue_ep6 = split (/\./, $sixthvalue_ep);
#my @splitvalue_ep7 = split (/\./, $seventhvalue_ep);
#
#my $firstdepth_ep = $splitvalue_ep1[0];
#my $seconddepth_ep = $splitvalue_ep2[0];
#my $thirddepth_ep = $splitvalue_ep3[0];
#my $fourthdepth_ep = $splitvalue_ep4[0];
#my $fifthdepth_ep = $splitvalue_ep5[0];
#my $sixthdepth_ep = $splitvalue_ep6[0];
#my $seventhdepth_ep = $splitvalue_ep7[0];
#
#$firstdepth_ep =~ s/\s//;
#$seconddepth_ep =~ s/\s//;
#$thirddepth_ep =~ s/\s//;
#$fourthdepth_ep =~ s/\s//;
#$fifthdepth_ep =~ s/\s//;
#$sixthdepth_ep =~ s/\s//;
#$seventhdepth_ep =~ s/\s//;
#
#my $min_ep = min grep {$_ > 300} $firstdepth_ep, $seconddepth_ep, $thirddepth_ep, $fourthdepth_ep, $fifthdepth_ep, $sixthdepth_ep, $seventhdepth_ep;
#
#if (defined($min_ep) && $min_ep >= 300){
#  $rardepth_ep = $min_ep;
#}else{
#  $rardepth_ep = 300;
#}
#
#close SUMEPI;
#
#my $rardepth_en;
#open (SUMENDO, $EndoSummary) or die "Can't open the summary file Endo, see line 388"."\n";
#my @data_en = <SUMENDO>;
#chomp @data_en;
#my $firstline_en = @data_en[16];
#my $secondline_en = @data_en[17];
#my $thirdline_en = @data_en[18];
#my $fourthline_en = @data_en[19];
#my $fifthline_en = @data_en[20];
#my $sixthline_en = @data_en[21];
#my $seventhline_en = @data_en[22];
#my @splitfirstline_en = split(/\:/,$firstline_en);
#my @splitsecondline_en = split(/\:/,$secondline_en);
#my @splitthirdline_en = split(/\:/,$thirdline_en);
#my @splitfourthline_en = split(/\:/,$fourthline_en);
#my @splitfifthline_en = split(/\:/,$fifthline_en);
#my @splitsixthline_en = split(/\:/,$sixthline_en);
#my @splitseventhline_en = split(/\:/,$seventhline_en);
#
#my $firstvalue_en = $splitfirstline_en[1];
#my $secondvalue_en = $splitsecondline_en[1];
#my $thirdvalue_en = $splitthirdline_en[1];
#my $fourthvalue_en = $splitfourthline_en[1];
#my $fifthvalue_en = $splitfifthline_en[1];
#my $sixthvalue_en = $splitsixthline_en[1];
#my $seventhvalue_en = $splitseventhline_en[1];
#
#my @splitvalue_en1 = split (/\./, $firstvalue_en);
#my @splitvalue_en2 = split (/\./, $secondvalue_en);
#my @splitvalue_en3 = split (/\./, $thirdvalue_en);
#my @splitvalue_en4 = split (/\./, $fourthvalue_en);
#my @splitvalue_en5 = split (/\./, $fifthvalue_en);
#my @splitvalue_en6 = split (/\./, $sixthvalue_en);
#my @splitvalue_en7 = split (/\./, $seventhvalue_en);
#
#my $firstdepth_en = $splitvalue_en1[0];
#my $seconddepth_en = $splitvalue_en2[0];
#my $thirddepth_en = $splitvalue_en3[0];
#my $fourthdepth_en = $splitvalue_en4[0];
#my $fifthdepth_en = $splitvalue_en5[0];
#my $sixthdepth_en = $splitvalue_en6[0];
#my $seventhdepth_en = $splitvalue_en7[0];
#
#$firstdepth_en =~ s/\s//;
#$seconddepth_en =~ s/\s//;
#$thirddepth_en =~ s/\s//;
#$fourthdepth_en =~ s/\s//;
#$fifthdepth_en =~ s/\s//;
#$sixthdepth_en =~ s/\s//;
#$seventhdepth_en =~ s/\s//;
#
#my $min_en = min grep {$_ > 500} $firstdepth_en, $seconddepth_en, $thirddepth_en, $fourthdepth_en, $fifthdepth_en, $sixthdepth_en, $seventhdepth_en;
#
#if (defined($min_en) && $min_en >= 500){
#  $rardepth_en = $min_en;
#}else{
#  $rardepth_en = 500;
#}
#close SUMENDO;
#
#sleep 2;
#
##future outputs:
#
#my $Rartable_Endo = $Rardir.$primer."_NO_0_taxa_cleaned_Endo_".$rardepth_en.".biom";
#my $Rartable_Epi = $Rardir.$primer."_NO_0_taxa_cleaned_Epi_".$rardepth_ep.".biom";
#my $Rartable_Wl = $Rardir.$primer."_NO_0_taxa_cleaned_Wl_".$rardepth_w.".biom";
#
#my $Rartable_Endo_t = $Rardir.$primer."_NO_0_taxa_cleaned_Endo_".$rardepth_en.".txt";
#my $Rartable_Epi_t = $Rardir.$primer."_NO_0_taxa_cleaned_Epi_".$rardepth_ep.".txt";
#my $Rartable_Wl_t = $Rardir.$primer."_NO_0_taxa_cleaned_Wl_".$rardepth_w.".txt";
#
####Now that we found the threshold for the rarefaction, and we purged out from the unwanted taxa, we rarefy accordingly.
#my $rarout1= $outdir."Rarefaction_Endo.out";
#my $rarout2= $outdir."Rarefaction_Wl.out";
#my $rarout3= $outdir."Rarefaction_Epi.out";
#
#my $rarcmd1 = "bsub -o $rarout1 single_rarefaction.py -i $TaxacleanEndo_b -o $Rartable_Endo -d $rardepth_en";
#my $rarcmd2 = "bsub -o $rarout2 single_rarefaction.py -i $TaxacleanWl_b -o $Rartable_Wl -d $rardepth_w";
#my $rarcmd3 = "bsub -o $rarout3 single_rarefaction.py -i $TaxacleanEpi_b -o $Rartable_Epi -d $rardepth_ep";
#
#my $run6 = system ('bash','-c',"$bsub && $qiime && $rarcmd1 && $rarcmd2 && $rarcmd3") == 0
#  or die "Can't perform the single rarefaction on the taxaclean table, see line 469"."\n";
#
#sleep 1 until -e $rarout1;
#sleep 1 until -e $rarout2;
#sleep 1 until -e $rarout3;
#
#
##now let's summarize by taxa, but let's keep in mind the differential formatting of Protists table and others:
#my $family;
#my $genus;
#my $species;
#
#if ($primer eq "PV4" || $primer eq "PV9"){
#   $family = 4;
#   $genus = 5;
#   $species = 6
#}else{
#   $family = 5;
#   $genus = 6;
#   $species = 7;
#}
#my $sumtaxout1= $outdir."TaxaSummary_Endo.out";
#my $sumtaxout2= $outdir."TaxaSummary_Wl.out";
#my $sumtaxout3= $outdir."TaxaSummary_Epi.out";
#
#my $summarycmd1= "bsub -o $sumtaxout1 summarize_taxa.py -i $Rartable_Endo -o $Taxasumdir -L $family,$genus,$species";
#my $summarycmd2= "bsub -o $sumtaxout2 summarize_taxa.py -i $Rartable_Wl -o $Taxasumdir -L $family,$genus,$species";
#my $summarycmd3= "bsub -o $sumtaxout3 summarize_taxa.py -i $Rartable_Epi -o $Taxasumdir -L $family,$genus,$species";
#
#my $run7 = system('bash', '-c', "$bsub && $qiime && $summarycmd1 && $summarycmd2 && $summarycmd3") == 0
#  or die "Can't get the summary (Endo, Wl, Epi) see line 499"."\n";
#
##Now we convert the cleaned otu tables, not summarized, in txt files for other usages, like loading in R.
#my $sndtxtconv1= $outdir."Txt_conversion_taxaclean_Endo.out";
#my $sndtxtconv2= $outdir."Txt_conversion_taxaclean_Wl.out";
#my $sndtxtconv3= $outdir."Txt_conversion_taxaclean_Epi.out";
#
#my $txtconvert1 = "bsub -o $sndtxtconv1 biom convert -b -i $TaxacleanEndo_b -o $TaxacleanEndo_t --table-type \"otu table\" --header-key taxonomy";
#my $txtconvert2 = "bsub -o $sndtxtconv2 biom convert -b -i $TaxacleanWl_b -o $TaxacleanWl_t --table-type \"otu table\" --header-key taxonomy";
#my $txtconvert3 = "bsub -o $sndtxtconv3 biom convert -b -i $TaxacleanEpi_b -o $TaxacleanEpi_t --table-type \"otu table\" --header-key taxonomy";
#
##The same with the rarefied ones, just in case
#my $sndtxtconv_rar1= $outdir."Txt_conversion_taxaclean_rar_Endo.out";
#my $sndtxtconv_rar2= $outdir."Txt_conversion_taxaclean_rar_Wl.out";
#my $sndtxtconv_rar3= $outdir."Txt_conversion_taxaclean_rar_Epi.out";
#
#my $txtconvert_rar1 = "bsub -o $sndtxtconv_rar1 biom convert -b -i $Rartable_Endo -o $Rartable_Endo_t --table-type \"otu table\" --header-key taxonomy";
#my $txtconvert_rar2 = "bsub -o $sndtxtconv_rar2 biom convert -b -i $Rartable_Wl -o $Rartable_Wl_t --table-type \"otu table\" --header-key taxonomy";
#my $txtconvert_rar3 = "bsub -o $sndtxtconv_rar3 biom convert -b -i $Rartable_Epi -o $Rartable_Epi_t --table-type \"otu table\" --header-key taxonomy";
#
#my $run8 = system('bash', '-c', "$bsub && $qiime && $txtconvert1 && $txtconvert2 && $txtconvert3 && $txtconvert_rar1 && $txtconvert_rar2 && $txtconvert_rar3") == 0
#  or die "Can't convert the final biom files into txt, see line 519"."\n";
#
#sleep 1 until -e $sndtxtconv1;
#sleep 1 until -e $sndtxtconv2;
#sleep 1 until -e $sndtxtconv3;
#sleep 1 until -e $sndtxtconv_rar1;
#sleep 1 until -e $sndtxtconv_rar2;
#sleep 1 until -e $sndtxtconv_rar3;

#Now create a merged OTU table Epi and Endo, in order to be used for plots in R such as alpha diversity, pcoa and so on
#outputs
#my $merged_table_b = $Rardir.$primer."_NO_0_taxa_cleaned_Epi_r_".$rardepth_ep."_Endo_r_".$rardepth_en."_merged.biom";
#my $merged_table_t = $Rardir.$primer."_NO_0_taxa_cleaned_Epi_r_".$rardepth_ep."_Endo_r_".$rardepth_en."_merged.txt";
#my $merge_out = $outdir."Merging_Epi_Endo.out";
#my $merge_out_t = $outdir."Conversion_in_Txt_merged.out";
#
#my $mergecmd = "bsub -o $merge_out merge_otu_tables.py -i $Rartable_Endo,$Rartable_Epi -o $merged_table_b";
#
#my $run9 = system('bash', '-c', "$bsub && $qiime && $mergecmd") == 0
#  or die "Can't create the merged Epi/Endo table, see line 539"."\n";
#
#sleep 1 until -e $merge_out;
#
#my $conv_merge ="bsub -o $merge_out_t biom convert -b -i $merged_table_b -o $merged_table_t --table-type \"otu table\" --header-key taxonomy";
#my $run10 = system('bash', '-c', "$bsub && $qiime && $conv_merge") == 0
#  or die "Can't convert the merged Epi/Endo table into txt, see line 545"."\n";
#
#sleep 1 until -e $merge_out_t;
#
##Let' create the summarised version of the merged
#
#my $mergesum_out = $outdir."Merged_table_summarising_taxa.out";
#my $summerge = "bsub -o $mergesum_out summarize_taxa.py -i $merged_table_b -o $Taxasumdir -L $family,$genus,$species";
#
#my $run11 = system('bash', '-c', "$bsub && $qiime && $summerge") == 0
#  or die "Can't taxa summarize the merged Epi/Endo table, see line 555"."\n";
#
#sleep 1 until -e $mergesum_out;

exit;
