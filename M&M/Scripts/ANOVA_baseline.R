#adonis sam
library(vegan)
#stats
#setwd("~/Desktop/Microbiome/final_ado/taxa_status/anova/")
args <- commandArgs(trailingOnly = TRUE)
mapfile <- "~/Desktop/Microbiome/Repository/M&M/Mapfile_rare_taxa.txt"
table <- "~/Desktop/Microbiome/Kingdom_effect/working_tables/BV5_merged_chao1.txt"#args[2]
#name <- args[2]
map1 <- read.csv(mapfile, header=T, row.names=NULL, sep = "\t", dec=".", stringsAsFactors = F)

div <- read.csv(table, row.names=NULL, header=T, dec=".", sep = "\t", stringsAsFactors = T)
divi <- data.frame(div, row.names=1)
head(divi)
divi$chao1 <- as.numeric(divi[,1])
#transform 
log_OTU <- log10(divi + 1)
log_OTU 
dim(log_OTU)
#we take only the samples present in the otu table
map2 <- map1#[,c(1:8)]
ordered_Map <- map2[match(rownames(log_OTU),map2$X.SampleID),]
head(ordered_Map)
dim(ordered_Map)
re_ordered_Map <- ordered_Map[complete.cases(ordered_Map),]
dim(log_OTU)
dim(re_ordered_Map)
rn <- rownames(log_OTU)[match(re_ordered_Map$X.SampleID,rownames(log_OTU))]
rn
re_sized_OTU <- data.frame(chao1=log_OTU[match(re_ordered_Map$X.SampleID,rownames(log_OTU)),])
head(re_sized_OTU)
rownames(re_sized_OTU) <- rn
dim(re_sized_OTU)
dim(re_ordered_Map)
head(re_sized_OTU)
head(re_ordered_Map)
dim(re_ordered_Map)
table(re_ordered_Map$Country)

re_ordered_Map$X.SampleID
rownames(re_sized_OTU)
re_ordered_Map$chao1 <- re_sized_OTU$chao1
#modify here
aov <- aov(chao1~ ***, re_ordered_Map) # insert here the variables
summary(aov)
ano <- anova(aov)
ano
aov

sink(paste("pos_anova.txt", sep = ""))
print(ano)
sink()

