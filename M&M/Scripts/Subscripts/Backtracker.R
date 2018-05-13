disclaimer <- "Hi, here you have your tracking of the network correlations, briefly for you to understand superquickly: the distance score is a transformed from pearson moment coefficient, 
              it means the same, is just on a different scale. 0 means superpositive correlation, 1 no correlation at all, and the closer you go to 2, this means negative interaction, important is then 
the distance class, all the distances were originally divided in classes in the mantel correlogram, here i kept the same rule (Sturge's rule), for making classes, this way, the distance class 1 
will be the same that you can find as the first dot in the correlogram and so on.
In other words, you should consider this file as a further explanation of the mantel correlogram, giving you more insight on WHO is making the correlation significant or not,
in order to do that, i just copied the pvalue from the correlogram in the pval column, but ATTENTION: that pvalue you should consider belonging to the CLASS, not to the single interaction"

args <- commandArgs(trailingOnly = TRUE)
mat1name <- args[1] 
mat1 <- args[2]
mat2name <- args[3]
mat2 <- args[4]
corr <- args[5]


mat <- read.csv(mat1, sep = "\t", row.names =1, header = TRUE)
mat_1 <- read.csv(mat2, sep = "\t", row.names =1, header = TRUE)
mat <- as.matrix(mat)
mat_1 <- as.matrix(mat_1)

map <- read.csv(corr, sep = "\t", row.names = NULL, header= TRUE, skip=7)
classes <- as.vector(map$Class.index)
n <- sum(map$Number.of.distances)
pval <- as.vector(map$p.value..Bonferroni.corrected.)

inc <- classes[3]-classes[2]
hinc <- inc/2
incr <- classes+hinc
decr <- classes-hinc
l <- length (classes)

for (i in 1:l){
if (i == 1){
fc <- which((mat_1 <incr[i] & mat_1>0), arr.ind = TRUE)  
}else{
fc<- which((mat_1 <incr[i] & mat_1>decr[i]), arr.ind = TRUE)
}
if (all(is.na(fc))==TRUE){
rown <- "No interactors in this class"
coln <- "No interactors in this class"
dint <- "NA"
counter <- 0
bonf <- "NA"
}else{
fc <- as.data.frame(fc)
names <- row.names(fc)
row <- fc[,1]
col <- fc[,2]
int <- mat_1[row,col]
dint <- diag(int)
rown <- row.names(int)
coln <- colnames(int)
counter <- c(1:length(rown))
bonf <- pval[i]
}
trackback <- data.frame(distance_ID = counter, interactor1 = rown, distance_score = dint, interactor2 = coln, distance_class= i, p_value = bonf)
write.table(trackback, append = TRUE, file= paste("Interaction_description_",mat1name,"_vs_",mat2name,".txt",sep = ""), sep= "\t", quote =F, row.names = F, dec =",")
#write.table(disclaimer, file= \"Readme.txt\", quote=F, row.names = F)
}