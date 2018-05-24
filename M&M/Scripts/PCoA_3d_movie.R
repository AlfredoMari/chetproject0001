library (vegan)
library(devtools)
library(vegan3d)
library(rgl)
library(magick)
args <- commandArgs(trailingOnly = TRUE)
f1 <- args[1] #feed mapfile
f2 <- args[2] #feed OTU table, no summarized
map1 <- read.csv(f1, header=T, row.names=NULL, sep = "\t", dec=".", stringsAsFactors = F)
head(map1)
div <- read.csv(f2, row.names=1,skip = 1, header=T, dec=".", sep = "\t", stringsAsFactors = F)
head(div)
no_tax <- function(x){
  x = x[,-which(colnames(x)=="taxonomy")]
}
divntx <- no_tax(div)
head(divntx)
map2 <- map1[match(colnames(divntx),map1$X.SampleID),]
tdiv<- t(divntx)
dim(tdiv)
dim(map2)

cap <- capscale(tdiv~Compartment+Country+Biome, map2, dist="bray")
scores(cap)
ordirgl(cap, display="species", choices= 1:3,type = "p",col= "green",arr.col="cadetblue1",radius=0.02, scaling= "symmetric")#c("blue", "green","black","magenta", "white","bisque","darkgoldenrod1","cadetblue1"), radius=0.03)
orglpoints(cap, display= "cn", col= c("black"), radius=0.02)
orglpoints(cap, display= "sites", col= c("deepskyblue"), radius=0.01)

dir <- "./pcoa/"
dir.create(dir)
play3d(spin3d(axis=c(-2,1,-1), rpm=3))
movie3d( spin3d(axis=c(-2,1,-1), rpm=1.5), duration = 40, fps= 15, movie= "protists_pcoa_40_15_minus_minus_final",type= "gif" ,
         dir = dir,convert=TRUE, clean=TRUE)
