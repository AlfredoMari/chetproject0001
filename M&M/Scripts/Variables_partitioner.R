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


#setwd("~/Desktop/Microbiome/Kingdom_effect/Lichens/Time_course/")     #provisional directory (to be fixed)
#we start with a general otu table, comprehensive of every marker clusater by super-sample ID
tab <- read.csv("~/Desktop/Microbiome/Kingdom_effect/working_tables/Time_course_data/Ftrad_otu_table_mc2_s2n50_L6.txt", sep = "\t", dec =".", skip =1, header = T, row.names=1)
map <- read.csv("~/Desktop/Microbiome/Kingdom_effect/working_tables/Time_course_data/mapfile_all_rn.txt", sep = "\t", dec =".", header = T, row.names=NULL)
#for sam's data we need to change the headers otherwise they don't match
colnames(tab) <- map$X.SampleID[match(colnames(tab),map$Barcode)]
colnames(tab)[ncol(tab)] <- "taxonomy" #re-insert the colnames taxonomy otherwise lost
head(tab)
naind <- which(complete.cases(colnames(tab))==TRUE)
tab <- tab[,c(naind)]
###
# divisore <- match(colnames(ctab),map$X.SampleID)
# d <-map$X.SampleID[divisore]
# comp_ser <- data.frame(samples = d, Compartment= map$Compartment[divisore])
# site_ser <- data.frame(samples = d, Sites = map$Site[divisore])

args  <- commandArgs( trailingOnly = TRUE )
kquery <- args[1]
#taxfile <- "~/Desktop/Microbiome/Kingdom_effect/BV5_OTU_table_merged_Endo4415_Epi504_L2.txt"
#div_file <- "~/Desktop/Microbiome/Kingdom_effect/BV5_merged_chao1.txt"
#e <- unlist(strsplit(taxfile, "_"))
name <- "BV5"

dir <- paste("./",kquery,"_on_",name,"/", sep= "")
#dir.create(dir)
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
# tab1 <- read.csv(taxfile, dec=".", sep = "\t", header=T,skip =1, row.names=1)
# taxa_tab <- t(tab1)
# indice <- match(rownames(taxa_tab), rownames(dfgr))
# New_tab_cue <- data.frame(taxa_tab, Cue_Status =dfgr$Cue_Status[indice], Site = dfgr$Sites[indice], Compartment= dfgr$Compartment[indice], stringsAsFactors=F)
# tax_sum_tot <- New_tab_cue[complete.cases(New_tab_cue),]
# tax_sum_epi <- tax_sum_tot[which(tax_sum_tot$Compartment=="R_Epi"),]
# tax_sum_endo <- tax_sum_tot[which(tax_sum_tot$Compartment=="R_Endo"),]
# 
# #indox <- match(rownames(New_tab_chloro), rownames(New_tab_ochro))
# #tax_sum_tot <- data.frame(New_tab_chloro, Red_Algae_Status =New_tab_ochro$Red_Algae_Status[indox], stringsAsFactors=F)
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# GetPalette <- colorRampPalette(cbbPalette, bias =1, interpolate = "spline", alpha = TRUE)
# #Epi_Site
# pos <- which(grepl("positive",tax_sum_epi$Cue_Status)==TRUE)
# neg <- which(grepl("negative",tax_sum_epi$Cue_Status)==TRUE)
# 
# freq_Site_pos <- table(tax_sum_epi$Site[pos])
# freq_Site_neg <- table(tax_sum_epi$Site[neg])
# freq_all <- rbind(freq_Site_pos, freq_Site_neg)
# #tu <- ddply (tax_sum_epi, ~ Site,numcolwise(sum))
# #ta <- data.frame(row.names=tu[,1],tu[,-1])
# tao <- apply(freq_all,1, function(x) x/sum(x))
# colnames(tao) <- c((paste(kquery,"_positive", sep = "")), (paste(kquery,"_negative", sep = "")))
# #te <- data.frame(ta,chloro_tax_sum$Green_Al)
# datm <- melt(tao)
# 
# colnames(datm) <- c("Site","Status","Frequence")
# ColorCount <- length(unique(datm$Site))
# 
# plots <- ggplot(datm, aes(x = Status, y = Frequence, fill = Site)) +
#   geom_bar(position = "stack", stat = "identity")+
#   #scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
#   scale_fill_manual(values = GetPalette(ColorCount))+
#   xlab("Grouping") + 
#   ylab("Abundance")+ 
#   theme (
#     panel.background = element_blank(),
#     legend.key = element_rect (fill = "white"),
#     legend.text = element_text (size = 20),
#     legend.title = element_text (size = 25, face = "bold"),
#     axis.text = element_text(size =20),
#     axis.line = element_line(color= "black", size =0.6),
#     axis.title = element_text(size =22, face = "bold") 
#     
#   )
# plotname01 <- "Final_plot_Compositional_Site_Epi3.pdf"
# ggsave (plots, file = paste(dir,plotname01,sep=""), width=20, height=9, units = "in", limitsize=FALSE)
# 
# #Endo_Site
# pos <- which(grepl("positive",tax_sum_endo$Cue_Status)==TRUE)
# neg <- which(grepl("negative",tax_sum_endo$Cue_Status)==TRUE)
# freq_Site_pos <- table(tax_sum_endo$Site[pos])
# freq_Site_neg <- table(tax_sum_endo$Site[neg])
# freq_all <- rbind(freq_Site_pos, freq_Site_neg)
# #tu <- ddply (tax_sum_epi, ~ Site,numcolwise(sum))
# #ta <- data.frame(row.names=tu[,1],tu[,-1])
# tao <- apply(freq_all,1, function(x) x/sum(x))
# colnames(tao) <- c((paste(kquery,"_positive", sep = "")), (paste(kquery,"_negative", sep = "")))
# #te <- data.frame(ta,chloro_tax_sum$Green_Al)
# datm <- melt(tao)
# colnames(datm) <- c("Site","Status","Frequence")
# ColorCount <- length(unique(datm$Site))
# 
# plots1 <- ggplot(datm, aes(x = Status, y = Frequence, fill = Site)) +
#   geom_bar(position = "stack", stat = "identity")+
#   #scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
#   scale_fill_manual(values = GetPalette(ColorCount))+
#   xlab("Grouping") + 
#   ylab("Abundance")+ 
#   theme (
#     panel.background = element_blank(),
#     legend.key = element_rect (fill = "white"),
#     legend.text = element_text (size = 20),
#     legend.title = element_text (size = 25, face = "bold"),
#     axis.text = element_text(size =20),
#     axis.line = element_line(color= "black", size =0.6),
#     axis.title = element_text(size =22, face = "bold") 
#     
#   )
# plotname01 <- "Final_plot_Compositional_Site_Endo3.pdf"
# ggsave (plots1, file = paste(dir,plotname01,sep=""), width=20, height=9, units = "in", limitsize=FALSE)
# 
# 
# 
# #Epi
# tu <- ddply (tax_sum_epi, ~ Cue_Status,numcolwise(sum))
# ta <- data.frame(row.names=tu[,1],tu[,-1])
# tao <- apply(ta,1, function(x) x/sum(x))
# #te <- data.frame(ta,chloro_tax_sum$Green_Al)
# datm <- melt(tao)
# 
# colnames(datm) <- c("Phylum","Status","Abundance")
# ColorCount <- length(unique(datm$Phylum))
# 
# plot1 <- ggplot(datm, aes(x = Status, y = Abundance, fill = Phylum)) +
#   geom_bar(position = "stack", stat = "identity")+
#   #scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
#   scale_fill_manual(values = GetPalette(ColorCount))+
#   xlab("Grouping") + 
#   ylab("Abundance")+ 
#   theme (
#     panel.background = element_blank(),
#     legend.key = element_rect (fill = "white"),
#     legend.text = element_text (size = 20),
#     legend.title = element_text (size = 25, face = "bold"),
#     axis.text = element_text(size =20),
#     axis.line = element_line(color= "black", size =0.6),
#     axis.title = element_text(size =22, face = "bold") 
#     
#   )
# plotname0 <- "Final_plot_Compositional_Epi3.pdf"
# ggsave (plot1, file = paste(dir,plotname0,sep=""), width=20, height=9, units = "in", limitsize=FALSE)
# 
# #Endo
# tu <- ddply (tax_sum_endo, ~ Cue_Status,numcolwise(sum))
# ta <- data.frame(row.names=tu[,1],tu[,-1])
# tao <- apply(ta,1, function(x) x/sum(x))
# #te <- data.frame(ta,chloro_tax_sum$Green_Al)
# datm <- melt(tao)
# 
# colnames(datm) <- c("Phylum","Status","Abundance")
# ColorCount <- length(unique(datm$Phylum))
# 
# plot2 <- ggplot(datm, aes(x = Status, y = Abundance, fill = Phylum)) +
#   geom_bar(position = "stack", stat = "identity")+
#   #scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
#   scale_fill_manual(values = GetPalette(ColorCount))+
#   xlab("Grouping") + 
#   ylab("Abundance")+ 
#   theme (
#     panel.background = element_blank(),
#     legend.key = element_rect (fill = "white"),
#     legend.text = element_text (size = 20),
#     legend.title = element_text (size = 25, face = "bold"),
#     axis.text = element_text(size =20),
#     axis.line = element_line(color= "black", size =0.6),
#     axis.title = element_text(size =22, face = "bold") 
#     
#   )
# plotname <- "Final_plot_Compositional_Endo3.pdf"
# ggsave (plot2, file = paste(dir,plotname,sep=""), width=20, height=9, units = "in", limitsize=FALSE)
# 
# #Alpha_diversity
# diver <- read.csv(div_file, sep = "\t", dec =".", header = T, row.names=1)  
# indice <- match(rownames(diver), rownames(dfgr))
# df <- data.frame(diver, Cue =dfgr$Cue_Status[indice], Compartment = dfgr$Compartment[indice], stringsAsFactors=F)
# df <- df[complete.cases(df),]
# df_epi <- data.frame(df[which(df$Compartment=="R_Epi"),])
# df_epi_to_print <- data.frame("Sample_ID" = rownames(df_epi),df_epi)
# write.table(df_epi_to_print, file = paste(dir,"Root_Epiphytes_alpha_div_source.txt",sep = ""), row.names= F, quote = F, dec =".", sep ="\t")
# df_endo <- data.frame(df[which(df$Compartment=="R_Endo"),])
# df_endo_to_print <- data.frame("Sample_ID" = rownames(df_endo),df_endo)
# write.table(df_endo_to_print, file = paste(dir,"Root_Endophytes_alpha_div_source.txt",sep = ""), row.names= F, quote = F, dec =".", sep ="\t")
# 
# #Epi
# df_epi<- data.frame(names = df_epi$Cue, chao1 =df_epi$chao1)
# colnames(df_epi) <- c(kquery, "chao1")
# datm <- melt(df_epi)
# colnames(datm) <- c("Sample","Metric","Alpha_Diversity")
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# GetPalette <- colorRampPalette(cbbPalette, bias =1, interpolate = "spline", alpha = TRUE)
# ColorCount <- length(unique(datm$Sample))
# 
# end <- ggplot(datm, aes(x = Sample, y = Alpha_Diversity))+#, fill = Metri)) +
#   geom_boxplot(lwd =1.5, fill = "white",aes(colour = Sample),outlier.shape=NA)+
#   geom_jitter(width=0.2,aes(colour = Sample, size=1))+
#   #scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
#   scale_colour_manual(values = GetPalette(ColorCount))+
#   xlab("Grouping") + 
#   ylab("Alpha diversity")+ 
#   theme (
#     panel.background = element_blank(),
#     legend.key = element_rect (fill = "white"),
#     legend.text = element_text (size = 20),
#     legend.title = element_text (size = 25, face = "bold"),
#     axis.text = element_text(size =20),
#     axis.line = element_line(color= "black", size =0.6),
#     axis.title = element_text(size =22, face = "bold") 
#     
#   )
# 
# plotname1 <- "Final_plot_alpha_Epi3.pdf"
# ggsave (end, file = paste(dir,plotname1,sep=""),device="pdf", width=20, height=9, units = "in")
# #now we calculate the stats
# histposep <- which(grepl("_positive",df_epi[,1]) == TRUE)
# histnegep <- which(grepl("_positive",df_epi[,1]) == FALSE)
# 
#   idx <- df_epi[,1]
#   levlen <- length(unique(idx))
#   if (levlen > 1){
#     wil <- wilcox.test(df_epi[,2]~df_epi[,1], data = df_epi, alternative = "two.sided")
#     wil$data.name <- paste("Bac_chao1 by ", kquery, sep = "")
#     graph_epi_neg <- df_epi$chao1[histnegep]
#     graph_epi_pos <- df_epi$chao1[histposep]
#   #print ("Yes")
#   }else{
#     wil <- "Cannot perform Wilcoxon/Mann-Withney test, there are in fact less than two classes to compare"
#     graph_epi_neg <- c(1)
#     graph_epi_pos <- c(1)
#   #print("No")
#     }
#   
# sink(paste(dir,"Wilcoxon_test_Epi.txt", sep = "")) 
# print(wil)
# sink()
# epi_distr_neg <- paste(dir,"Distribution_alpha_div_Bac_Epi_Neg.pdf", sep ="")
# pdf(epi_distr_neg)
# hist(graph_epi_neg)
# dev.off()
# epi_distr_pos <- paste(dir,"Distribution_alpha_div_Bac_Epi_Pos.pdf", sep ="")
# pdf(epi_distr_pos)
# hist(graph_epi_pos)
# dev.off()
# 
# 
# 
# #Endo
# df_endo<- data.frame(names = df_endo$Cue, chao1 =df_endo$chao1)
# colnames(df_endo) <- c(kquery, "chao1")
# datm <- melt(df_endo)
# colnames(datm) <- c("Sample","Metric","Alpha_Diversity")
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# GetPalette <- colorRampPalette(cbbPalette, bias =1, interpolate = "spline", alpha = TRUE)
# ColorCount <- length(unique(datm$Sample))
# 
# end <- ggplot(datm, aes(x = Sample, y = Alpha_Diversity))+#, fill = Metri)) +
#   geom_boxplot(lwd =1.5, fill = "white",aes(colour = Sample),outlier.shape=NA)+
#   geom_jitter(width=0.2,aes(colour = Sample, size=1))+
#   #scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
#   scale_colour_manual(values = GetPalette(ColorCount))+
#   xlab("Grouping") + 
#   ylab("Alpha diversity")+ 
#   theme (
#     panel.background = element_blank(),
#     legend.key = element_rect (fill = "white"),
#     legend.text = element_text (size = 20),
#     legend.title = element_text (size = 25, face = "bold"),
#     axis.text = element_text(size =20),
#     axis.line = element_line(color= "black", size =0.6),
#     axis.title = element_text(size =22, face = "bold") 
#     
#   )
# 
# plotname2 <- "Final_plot_alpha_Endo.pdf"
# ggsave (end, file = paste(dir,plotname2,sep=""),device="pdf", width=20, height=9, units = "in")
# #now we calculate the stats
# histpos <- which(grepl("_positive",df_endo[,1]) == TRUE)
# histneg <- which(grepl("_positive",df_endo[,1]) == FALSE)
# 
# idx1 <- df_endo[,1]
# levlen1 <- length(unique(idx1))
#   if (levlen1 > 1){
#     wil1 <- wilcox.test(df_endo[,2]~df_endo[,1], data = df_endo, alternative = "two.sided")
#     wil1$data.name <- paste("Bac_chao1 by ", kquery, sep = "")
#     graph_endo_neg <- df_endo$chao1[histneg]
#     graph_endo_pos <- df_endo$chao1[histpos]
#   }else{
#     wil1 <- "Cannot perform Wilcoxon/Mann-Withney test, there are in fact less than two classes to compare"
#     graph_endo_neg <- c(1)
#     graph_endo_pos <- c(1)
#   }
# 
# sink(paste(dir,"Wilcoxon_test_Endo.txt", sep = ""))   
# print(wil1)
# sink()
# 
# endo_distr_neg <- paste(dir,"Distribution_alpha_div_Bac_Endo_Neg.pdf", sep ="")
# pdf(endo_distr_neg)
# hist(graph_endo_neg)
# dev.off()
# 
# endo_distr_pos <- paste(dir,"Distribution_alpha_div_Bac_Endo_Pos.pdf", sep ="")
# pdf(endo_distr_pos)
# hist(graph_endo_pos)
# dev.off()
# 
# map
# ep <- which(grepl("R_",map$Compartment)==TRUE)
# root_map <- map[ep,]
# write.table(root_map, file = "Root_mapfile.txt", sep = "\t", quote = F, row.names=F, dec= ".")
