#setwd("/Volumes/biodata/dep_psl/grp_kemen/working_directory_alfredo/New_Methods/lobo_compar/lobo_250-8/")
args <- commandArgs(trailingOnly = T)
net1 <- args[1]
net2 <- args[2]
n1 <- read.csv(net1, sep = ",", header=T, row.names=NULL, dec=".", stringsAsFactors = F)
n2 <- read.csv(net2, sep = ",", header=T, row.names=NULL, dec=".", stringsAsFactors = F)

head(n1)
head(n2)
listn1 <- strsplit(as.character(n1$shared.name), "->")
listn2 <- strsplit(as.character(n2$shared.name), "->")


or1 <- unlist(listn1)[c(TRUE,FALSE)]
tgt1 <- unlist(listn1)[c(FALSE,TRUE)]
n1$origin_interactor <- or1
n1$target_interactor <- tgt1
write.table(n1,file ="minus_edge_table_div.txt", sep ="\t", quote =F, row.names=F, col.names=T)

or2 <- unlist(listn2)[c(TRUE,FALSE)]
tgt2 <- unlist(listn2)[c(FALSE,TRUE)]
n2$origin_interactor <- or2
n2$target_interactor <- tgt2
write.table(n2,file ="plus_edge_table_div.txt", sep ="\t", quote =F, row.names=F, col.names=T)
#head(n2)


