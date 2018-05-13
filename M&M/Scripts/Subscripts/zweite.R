args <- commandArgs(trailingOnly = TRUE)
args[1] <- Itrack1 
args[2] <- Itrack2 

track1 <- read.csv(Itrack1, sep ="\t", header = F, row.names = NULL)
track2 <- read.csv(Itrack2, sep ="\t", header = F, row.names = NULL)
                    ID <- as.vector(track1[,1])
                    Int <- as.vector(track1[,2])
                    dist <- as.vector(track1[,3])
                    Int2 <- as.vector(track1[,4])
                    distclass <- as.vector(track1[,5])
                    pval <- as.vector(track1[,6])
                    
                    sign <- which(pval<0.05)
                    
                    significant <- data.frame(distance_ID = ID[sign],Interactor1 =Int[sign],distance_score= dist[sign],Interactor2 =Int2[sign],Distance_class =distclass[sign],P_value=pval[sign], Original_matrix = 2)
                    
                    ID2 <- as.vector(track2[,1])
                    Int2 <- as.vector(track2[,2])
                    dist2 <- as.vector(track2[,3])
                    Int22 <- as.vector(track2[,4])
                    distclass2 <- as.vector(track2[,5])
                    pval2 <- as.vector(track2[,6])
                    
                    sign2 <- which(pval2<0.05)
                    
                    significant2 <- data.frame(distance_ID = ID2[sign2],Interactor1 =Int2[sign2],distance_score= dist2[sign2],Interactor2 =Int22[sign2],Distance_class =distclass2[sign2],P_value=pval2[sign2], Original_matrix =1)
                    
                    merge <- rbind(significant, significant2)
                    
                    core_merged <- unique(merge[,2:4], incomparables = F)
                    
                    back_merge_Int <- as.vector(unique(match(significant[,2], core_merged[,1])))
                    back_merge_dist <- as.vector(unique(match(significant[,3], core_merged[,2])))
                    back_merge_Int2 <- as.vector(unique(match(significant[,4],core_merged[,3])))
                    
                    
                    coordinates <- c(back_merge_Int,back_merge_dist,back_merge_Int2)
                    u_coordinates <- unique(coordinates)
                    
                    back_merge_Int2 <- as.vector(unique(match(significant2[,2], core_merged[,1])))
                    back_merge_dist2 <- as.vector(unique(match(significant2[,3], core_merged[,2])))
                    back_merge_Int22<- as.vector(unique(match(significant2[,4],core_merged[,3])))
                    
                    
                    coordinates2 <- c(back_merge_Int2,back_merge_dist2,back_merge_Int22)
                    u_coordinates2 <- unique(coordinates2)
                    
                    final_core_1 <- as.data.frame(list(significant[u_coordinates,]))
                    final_core_2 <- as.data.frame(list(significant2[u_coordinates2,]))
                    final_core <- rbind(final_core_1,final_core_2)
                    
                    write.table(final_core, file = paste("Significant_Interactions_whole_both_net_merged_core.txt",sep = "\t"), col.names = T, row.names = F, quote = F, dec = ",", sep= "\t")
                    