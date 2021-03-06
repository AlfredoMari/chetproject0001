#Read in the first edge table
setwd("/Volumes/biodata/dep_psl/grp_kemen/working_directory_alfredo/New_Methods/brac2/")
args <- commandArgs(trailingOnly = T)
fstNetName <- args[1]
fstNet <- "minus_edge_table_div.txt"
sndNetName <- args[3]
sndNet <- args[4]
data  <- read.csv(fstNet, sep ="\t", header= T, dec =".")
                   #Read in the second edge table
                   data_1  <- read.csv(sndNet, sep ="\t", header= T, dec =".")
                   
                   #Extract the origin interactor, the target, and the rsquare value
                   org  <- data$origin_interactor
                   tgt  <- data$target_interactor
                   rsq  <- data$weight
                 #  head(format(rsq, scientific=F))
                   #rsqp <- log(rsq)/10
                   #new_rsq <- (2-2^(rsq))   #This will keep the original values in, just in a different scale, values from 0 to 1 indicate positive interaction, 1 is no interaction, from 1 to 2 is negative interaction
                   new_rsq <- rsq
                  # head(new_rsq)
                   #Feed all the raw values in a data.frame, this will be the source for matching the rsq values in the matrix
                   the_dataframe  <- data.frame(org,new_rsq,tgt)
                   #Remove the doubles for origin and target
                   correct_tgt  <- unique(tgt,incomparables= F)
                   correct_org  <- unique(org,incomparables= F)
                   
                   #the same for the second set
                   org_1  <- data_1$origin_interactor
                   tgt_1  <- data_1$target_interactor
                   rsq_1  <- data_1$weight
                   #rsqp_1 <- log(rsq_1)/10
                #   new_rsq_1 <- (2-2^(rsq_1))   
                   new_rsq_1 <- rsq_1
                   the_dataframe_1  <- data.frame(org_1,new_rsq_1,tgt_1)
                   correct_tgt_1  <- unique(tgt_1,incomparables= F)
                   correct_org_1  <- unique(org_1,incomparables= F)
                   
                   ###BUILDING THE OUTPUT MATRIX###
                   #Merge the origin and the interactor for set one
                   vec  <- c(as.character(correct_org), as.character(correct_tgt))
                   vec1  <- unique(vec,incomparables= F)
                   #the same with second set
                   vec_1  <- c(as.character(correct_org_1), as.character(correct_tgt_1))
                   vec1_1  <- unique(vec_1,incomparables= F)
                   #Merge all together
                   all_interactors  <- c(vec1, vec1_1)
                   all_interactors_end  <- unique(all_interactors, incomparables=F)
                   l  <- length(all_interactors_end)
                   #Output matrices, they are made by the same rownames and colnames, so the size is the same.
                   #1
                   dm = matrix( nrow= l, ncol=l, dimnames = list(all_interactors_end, all_interactors_end))
                   #2
                   dm_1 = matrix( nrow= l, ncol=l, dimnames = list(all_interactors_end, all_interactors_end))
                   
                   ###FILLING THE OUTPUT MATRIX###
                   #Goes back into the rawdata and for every interaction, it fills with the correspondant adj_rsq value
                   #1st Matrix
                   for (i in 1:nrow(the_dataframe)){ #iterate through lines of input file
                   
                   dm[match(the_dataframe[i,1],rownames(dm)), match(the_dataframe[i,3],colnames(dm))] = the_dataframe[i,2]
                   
                   }
                   dm[is.na(dm)] <- 1  #converts the NA into 1 (0 would mean superpositive interaction, this way we assign the NA to the non correlated, makes more sense and prevents a positive interaction inflation, that was the bug of the previous version)
                   
                   #2nd Matrix
                   for (j in 1:nrow(the_dataframe_1)){ #iterate through lines of input file
                   
                   dm_1[match(the_dataframe_1[j,1],rownames(dm_1)), match(the_dataframe_1[j,3],colnames(dm_1))] = the_dataframe_1[j,2]
                   
                   }
                   dm_1[is.na(dm_1)] <- 1 #converts the NA into 1
                   #NOTICE: the filling step takes time.
                   
                   ###PREPARING THE MANTEL TEST###
                   #converts the two matrices into distance matrices, so that we obtained a sparse and hollow matrix.
                   edm <- as.dist(dm)
                   edm_1 <- as.dist(dm_1)
                   
                   #Check the size, in order to continue, edm and edm_1 should have the same size
                   attr(edm, "Size")
   attr(edm_1,"Size")

#A dist. object cannot be printed, that is why we convert it into matrix before writing it down
   hol_edm <- as.matrix(edm)
   hol_edm_1  <- as.matrix(edm_1)

#writes down the distance matrices which can be feeded into compare_distance_matrices.py for mantel test
   write.table(hol_edm, file =paste(fstNetName,"_hollow_dist.txt",sep = ""), sep= "\t", col.names =NA, row.names=TRUE, quote =F)
   write.table(hol_edm_1, file =paste(sndNetName,"_hollow_dist.txt",sep = ""), sep= "\t", col.names =NA, row.names=TRUE, quote =F)

#Further utilities: write the original matrices down for testing them elsewhere:
   write.table(dm, file =paste(fstNetName,"_row_matrix.txt",sep=""), sep= "\t", col.names =NA, row.names=TRUE, quote =F)
   write.table(dm_1, file =paste(sndNetName,"_row_matrix.txt",sep = ""), sep= "\t", col.names =NA, row.names=TRUE, quote =F)                                                                          