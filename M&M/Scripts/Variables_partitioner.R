library ('gplots')
library ('vegan')
library ('ggplot2')
library ('RColorBrewer')
library ('MASS')
library ('reshape2')
library ('RColorBrewer')
library('gdata')
library('plyr')
library('dplyr')


#we start with a general otu table, comprehensive of every marker clusater by super-sample ID
args <- commandArgs(trailingOnly= TRUE)
f1 <- args[1] # introduce summarized taxa otu table only
f2 <- args[2] # introduce mapfile
tab <- read.csv(f1, sep = "\t", dec =".", skip =1, header = T, row.names=1) 
map <- read.csv(f2, sep = "\t", dec =".", header = T, row.names=NULL) 
#for common garden data we need to change the headers otherwise they don't match
colnames(tab) <- map$X.SampleID[match(colnames(tab),map$Barcode)]
colnames(tab)[ncol(tab)] <- "taxonomy" #re-insert the colnames taxonomy otherwise lost
head(tab)
naind <- which(complete.cases(colnames(tab))==TRUE)
tab <- tab[,c(naind)]

args  <- commandArgs( trailingOnly = TRUE )
kquery <- args[1]
name <- "BV5"

dir <- paste("./",kquery,"_on_",name,"/", sep= "")
dir.create(dir)
ctab <- tab[,-which(names(tab)=="taxonomy")]
cue <- grepl(kquery,rownames(tab))    #select from the table the otus assigned to cuerophyta

cue_index <- which(cue == "TRUE")

cue_sub <- tab[cue_index,]  #subset of only cuerophyta otus

cue_sum <- colSums(cue_sub[,-which(colnames(cue_sub)=="taxonomy")])  # sum all the abundances per sample 

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

#ind_comp <- match(rownames(dfgr),comp_ser$samples)
#ind_sites <- match(rownames(dfgr), site_ser$samples)
#dfgr <- data.frame(dfgr,Compartment =comp_ser$Compartment[ind_comp],Sites = site_ser$Sites[ind_sites], stringsAsFactors = F)
dfgr
New_mapfile <-data.frame(map, Taxa_status = dfgr$Cue_Status[match(map$X.SampleID,rownames(dfgr))])
#names(New_mapfile$Taxa_status) <- paste(kquery,"_staus", sep = "")
names(New_mapfile)[names(New_mapfile) == "Taxa_status"] <- paste(kquery,"_taxa_status", sep = "")
head(New_mapfile)
write.table(New_mapfile,file ="Mapfile_sam.txt", quote=F, dec= ".", sep = "\t", row.names= F, col.names = T)
