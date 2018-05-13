library(Hmisc) # cut2
library ('gplots')
library ('vegan')
library ('ggplot2')
library ('RColorBrewer')
library ('MASS')
library ('reshape2')
library ('RColorBrewer')
#install.packages('doBy')
library('gdata')
library('plyr')
library('dplyr')


#setwd("~/Desktop/Microbiome//Kingdom_effect")     #provisional directory (to be fixed)
#we start with a general otu table, comprehensive of every marker clusater by super-sample ID
tab <- read.csv("~/Desktop/Microbiome/Kingdom_effect/working_tables/PV4_OTU_table_merged_Endo3055_Epi412.txt", sep = "\t", dec =".", skip =1, header = T, row.names=1)
ctab <- tab[,-which(names(tab)=="taxonomy")]
map <- read.csv("~/Desktop/Master_mapfile_RPCOA_complete_4Net.txt", sep = "\t", dec =".", header = T, row.names=NULL)

divisore <- match(colnames(ctab),map$X.SampleID)
d <-map$X.SampleID[divisore]
comp_ser <- data.frame(samples = d, Compartment= map$Compartment[divisore])
site_ser <- data.frame(samples = d, Sites = map$Site[divisore])

#e <- unlist(strsplit(taxfile, "_"))
args  <- commandArgs( trailingOnly = TRUE )
kquery <- args[1]
div_choice <- args[2]
if(div_choice == "B"){
  div_file <- "~/Desktop/Microbiome/Kingdom_effect/working_tables/BV5_merged_chao1.txt"
}else if(div_choice =="F"){
  div_file <-   "~/Desktop/Microbiome/Kingdom_effect/working_tables/FITS2_merged_chao1.txt"
}else if(div_choice =="O"){
  div_file <- "~/Desktop/Microbiome/Kingdom_effect/working_tables/Otrad_merged_chao1.txt" 
}
name <- unlist(strsplit(unlist(strsplit(div_file,"\\/"))[6],"_"))[1]


dir <- paste("./",kquery,"_taxa_only_on_",name,"/", sep= "")
dir.create(dir)

cue <- grepl(kquery,tab$taxonomy)    #select from the table the otus assigned to cuerophyta

cue_index <- which(cue == "TRUE")

cue_sub <- tab[cue_index,-which(names(tab) == "taxonomy")]  #subset of only cuerophyta otus

cue_sum <- colSums(cue_sub)  # sum all the abundances per sample 

cue_pos <- which(cue_sum > 0, arr.ind = TRUE)  #identifies samples which have at least one otus which includes cuero
cue_neg <- which(cue_sum == 0, arr.ind = TRUE)  #identifies samples in which cuerophyta are completely absent
cue_pos
cue_neg
w_cue <- cue_sub[,cue_pos]  #from the cuero_subset, picks up only the samples with at least one otus assigned to cuero
wo_cue <- cue_sub[,cue_neg]  #from the cuero_subset, picks up only the samples in which we have no otus assigned to cuero

w_cuein <- colnames(w_cue) #takes the sample IDs from the positive cuero subset
wo_cuein <- colnames(wo_cue) #takes the sample IDs from the 0 cuero subset

#we add a description stating if each sample is positive/negative for red algae or green algae, adding the sum value of the whole cuerophyta abundance per sample
#under the column Abundance
wcue_df <- data.frame(ID= w_cuein, Cue_Status = paste(kquery,"_positive", sep = ""), Cue_Abundance = colSums(w_cue), row.names=1)
wocue_df <- data.frame(ID= wo_cuein, Cue_Status= paste(kquery,"_negative", sep = ""),Cue_Abundance = colSums(wo_cue), row.names=1)
dfgr <- rbind(wcue_df,wocue_df)

ind_comp <- match(rownames(dfgr),comp_ser$samples)
ind_sites <- match(rownames(dfgr), site_ser$samples)
dfgr <- data.frame(dfgr,Compartment =comp_ser$Compartment[ind_comp],Sites = site_ser$Sites[ind_sites], stringsAsFactors = F)
dfgr


#indox <- match(rownames(New_tab_chloro), rownames(New_tab_ochro))
#tax_sum_tot <- data.frame(New_tab_chloro, Red_Algae_Status =New_tab_ochro$Red_Algae_Status[indox], stringsAsFactors=F)


#Alpha_diversity
diver <- read.csv(div_file, sep = "\t", dec =".", header = T, row.names=1)  
indice <- match(rownames(diver), rownames(dfgr))
df <- data.frame(diver, Cue =dfgr$Cue_Status[indice], Compartment = dfgr$Compartment[indice], stringsAsFactors=F)
df <- df[complete.cases(df),]
df_epi <- data.frame(df[which(df$Compartment=="Epi"),])
df_epi_to_print <- data.frame("Sample_ID" = rownames(df_epi),df_epi)
write.table(df_epi_to_print, file = paste(dir,"Leaf_Epiphytes_alpha_div_source.txt",sep = ""), row.names= F, quote = F, dec =".", sep ="\t")
df_endo <- data.frame(df[which(df$Compartment=="Endo"),])
df_endo_to_print <- data.frame("Sample_ID" = rownames(df_endo),df_endo)
write.table(df_endo_to_print, file = paste(dir,"Leaf_Endophytes_alpha_div_source.txt",sep = ""), row.names= F, quote = F, dec =".", sep ="\t")

#Epi
df_epi<- data.frame(names = df_epi$Cue, chao1 =df_epi$chao1)
colnames(df_epi) <- c(kquery, "chao1")
datm <- melt(df_epi)
colnames(datm) <- c("Sample","Metric","Alpha_Diversity")
dati <- aggregate(datm$Alpha_Diversity, list(datm$Sample), mean)
if(nrow(dati)>1){
  scarto <- (dati[2,2]-dati[1,2])
  if(scarto >0){}else{scarto <- scarto*-1}
  #normalize it
  scar <- (scarto*100)/min(dati$x)
  dati1 <- t(dati)
  dat <- data.frame(dati, perc_gap = c("-",scar))
}else{
  dat <- "Can't compute a gap, there is infact only one group present"
}
write.table(dat, file= paste(dir,kquery,"_on_",name,"_fold_alpha_diversity_Epi.txt", sep = ""), quote = F, dec = ".", sep = "\t", row.names = F)

#now we calculate the stats
histposep <- which(grepl("_positive",df_epi[,1]) == TRUE)
histnegep <- which(grepl("_positive",df_epi[,1]) == FALSE)

  idx <- df_epi[,1]
  levlen <- length(unique(idx))
  if (levlen > 1){
    wil <- wilcox.test(df_epi[,2]~df_epi[,1], data = df_epi, alternative = "two.sided")
    wil$data.name <- paste("Bac_chao1 by ", kquery, sep = "")
    graph_epi_neg <- df_epi$chao1[histnegep]
    graph_epi_pos <- df_epi$chao1[histposep]
  #print ("Yes")
  }else{
    wil <- "Cannot perform Wilcoxon/Mann-Withney test, there are in fact less than two classes to compare"
    graph_epi_neg <- c(1)
    graph_epi_pos <- c(1)
  #print("No")
    }
  
sink(paste(dir,"Wilcoxon_test_Epi.txt", sep = "")) 
print(wil)
sink()
epi_distr_neg <- paste(dir,"Distribution_alpha_div_Bac_Epi_Neg.pdf", sep ="")
pdf(epi_distr_neg)
hist(graph_epi_neg)
dev.off()
epi_distr_pos <- paste(dir,"Distribution_alpha_div_Bac_Epi_Pos.pdf", sep ="")
pdf(epi_distr_pos)
hist(graph_epi_pos)
dev.off()



#Endo
df_endo<- data.frame(names = df_endo$Cue, chao1 =df_endo$chao1)
colnames(df_endo) <- c(kquery, "chao1")
datm <- melt(df_endo)
colnames(datm) <- c("Sample","Metric","Alpha_Diversity")
dati <- aggregate(datm$Alpha_Diversity, list(datm$Sample), mean)
if(nrow(dati)>1){
  scarto <- (dati[2,2]-dati[1,2])
  if(scarto >0){}else{scarto <- scarto*-1}
  #normalize it
  scar <- (scarto*100)/min(dati$x)
  dati1 <- t(dati)
  dat <- data.frame(dati, perc_gap = c("-",scar))
}else{
  dat <- "Can't compute a gap, there is infact only one group present"
}
write.table(dat, file= paste(dir,kquery,"_on_",name,"_fold_alpha_diversity_Endo.txt", sep = ""), quote = F, dec = ".", sep = "\t", row.names = F)


#now we calculate the stats
histpos <- which(grepl("_positive",df_endo[,1]) == TRUE)
histneg <- which(grepl("_positive",df_endo[,1]) == FALSE)

idx1 <- df_endo[,1]
levlen1 <- length(unique(idx1))
  if (levlen1 > 1){
    wil1 <- wilcox.test(df_endo[,2]~df_endo[,1], data = df_endo, alternative = "two.sided")
    wil1$data.name <- paste("Bac_chao1 by ", kquery, sep = "")
    graph_endo_neg <- df_endo$chao1[histneg]
    graph_endo_pos <- df_endo$chao1[histpos]
  }else{
    wil1 <- "Cannot perform Wilcoxon/Mann-Withney test, there are in fact less than two classes to compare"
    graph_endo_neg <- c(1)
    graph_endo_pos <- c(1)
  }

sink(paste(dir,"Wilcoxon_test_Endo.txt", sep = ""))   
print(wil1)
sink()

endo_distr_neg <- paste(dir,"Distribution_alpha_div_Bac_Endo_Neg.pdf", sep ="")
pdf(endo_distr_neg)
hist(graph_endo_neg)
dev.off()

endo_distr_pos <- paste(dir,"Distribution_alpha_div_Bac_Endo_Pos.pdf", sep ="")
pdf(endo_distr_pos)
hist(graph_endo_pos)
dev.off()

