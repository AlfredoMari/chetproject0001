#adonis sam
library(vegan)
#stats
setwd("~/Desktop/Microbiome/final_ado/taxa_status/")
args <- commandArgs(trailingOnly = TRUE)
mapfile <- "~/Desktop/Microbiome/Kingdom_effect/Mapfile_rare_taxa.txt"
table <- "~/Desktop/Microbiome/Kingdom_effect/working_tables/BV5_merged_chao1.txt"#args[2]
#name <- args[2]
map1 <- read.csv(mapfile, header=T, row.names=NULL, sep = "\t", dec=".", stringsAsFactors = F)
#map1 <- map1[,c(1:23)]
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
#modify here
ado2 <- adonis(re_sized_OTU~ Lobosa_taxa_status *      Bracteacoccus_taxa_status * Chlorophyta_taxa_status * Ciliophora_taxa_status *      Metazoa_taxa_status * Ochrophyta_taxa_status * Cercozoa_taxa_status*Discoba_taxa_status*Year, data = re_ordered_Map, permutations = 999, method = "bray", strata =NULL)
ado2
#and here
sink(paste("pos_taxa_status_bac.txt", sep = ""))
print(ado2)
sink()


