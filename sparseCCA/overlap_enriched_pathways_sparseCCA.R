## Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu

## This script identifies pathway overlaps across diseases.

## Initial setup 
rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
setwd("/path/to/working_directory")

################ Input case-specific for each disease ############

msigdb_collection_name <- "C2_CP/case_minus_control_kegg_pid/FDR_cutoff/25_300_5/"

## CRC
dataset <- "CRC/case"
filename <- paste0("CRC_case_specific_CompGO_FDR_0.1.txt")
filepath <- paste0("./data/analysis/",dataset,"/output_sparseCCA/grid_search/enrichment_sparseCCA_msigDB/",msigdb_collection_name,"/",filename)
CRC_case_minus_control_FDR_0.1 <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
all(CRC_case_minus_control_FDR_0.1$p_adj < 0.1)
dim(CRC_case_minus_control_FDR_0.1)


## IBD
dataset <- "IBD/case"
filename <- paste0("IBD_case_specific_CompGO_FDR_0.1.txt")
filepath <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/enrichment_sparseCCA_msigDB/",msigdb_collection_name,"/",filename)
IBD_case_minus_control_FDR_0.1 <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
dim(IBD_case_minus_control_FDR_0.1) 
all(IBD_case_minus_control_FDR_0.1$p_adj < 0.1)

## IBS
dataset <- "IBS/case"
filename <- paste0("IBS_case_specific_CompGO_FDR_0.1.txt")
filepath <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/enrichment_sparseCCA_msigDB/",msigdb_collection_name,"/",filename)
IBS_case_minus_control_FDR_0.1 <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
dim(IBS_case_minus_control_FDR_0.1) 

CRC_case_minus_control <- CRC_case_minus_control_FDR_0.1
IBD_case_minus_control <- IBD_case_minus_control_FDR_0.1
IBS_case_minus_control <- IBS_case_minus_control_FDR_0.1


################ Compute overlaps across diseases ###################

## (CRC ^ IBD ^ IBS)
CRC_IBD_IBS_case_minus_control_overlap <- Reduce(function(x,y) merge (x,y, by = "pathway" ), 
                                                 list(CRC_case_minus_control, IBD_case_minus_control, IBS_case_minus_control))
dim(CRC_IBD_IBS_case_minus_control_overlap) 

## sort by adjusted pval(CRC, IBD, IBS)
CRC_IBD_IBS_case_minus_control_overlap <- CRC_IBD_IBS_case_minus_control_overlap[order(CRC_IBD_IBS_case_minus_control_overlap$p_adj.x,
                                                                                       CRC_IBD_IBS_case_minus_control_overlap$p_adj.y,
                                                                                       CRC_IBD_IBS_case_minus_control_overlap$p_adj),]



## (CRC ^ IBD) = (CRC ^ IBD) - (CRC ^ IBD ^ IBS)
CRC_IBD_case_minus_control_overlap <- merge(CRC_case_minus_control, IBD_case_minus_control, by = "pathway" )
dim(CRC_IBD_case_minus_control_overlap) 
CRC_IBD_case_minus_control_overlap <- CRC_IBD_case_minus_control_overlap[!(CRC_IBD_case_minus_control_overlap$pathway %in% CRC_IBD_IBS_case_minus_control_overlap$pathway),]
dim(CRC_IBD_case_minus_control_overlap) 

## sort
CRC_IBD_case_minus_control_overlap <- CRC_IBD_case_minus_control_overlap[order(CRC_IBD_case_minus_control_overlap$p_adj.x,
                                                                               CRC_IBD_case_minus_control_overlap$p_adj.y),]


## (CRC ^ IBS) = (CRC ^ IBS) - (CRC ^ IBD ^ IBS)
CRC_IBS_case_minus_control_overlap <- merge(CRC_case_minus_control, IBS_case_minus_control, by = "pathway" )
dim(CRC_IBS_case_minus_control_overlap) 
CRC_IBS_case_minus_control_overlap <- CRC_IBS_case_minus_control_overlap[!(CRC_IBS_case_minus_control_overlap$pathway %in% CRC_IBD_IBS_case_minus_control_overlap$pathway),]
dim(CRC_IBS_case_minus_control_overlap) 

##sort
CRC_IBS_case_minus_control_overlap <- CRC_IBS_case_minus_control_overlap[order(CRC_IBS_case_minus_control_overlap$p_adj.x,
                                                                               CRC_IBS_case_minus_control_overlap$p_adj.y),]

## (IBD ^ IBS) = (IBD ^ IBS) - (IBD ^ IBD ^ IBS)
IBD_IBS_case_minus_control_overlap <- merge(IBD_case_minus_control, IBS_case_minus_control, by = "pathway" )
dim(IBD_IBS_case_minus_control_overlap) 
IBD_IBS_case_minus_control_overlap <- IBD_IBS_case_minus_control_overlap[!(IBD_IBS_case_minus_control_overlap$pathway %in% CRC_IBD_IBS_case_minus_control_overlap$pathway),]
dim(IBD_IBS_case_minus_control_overlap) 

##sort
IBD_IBS_case_minus_control_overlap <- IBD_IBS_case_minus_control_overlap[order(IBD_IBS_case_minus_control_overlap$p_adj.x,
                                                                               IBD_IBS_case_minus_control_overlap$p_adj.y),]



## CRC_specific = CRC - (CRC ^ IBD) - (CRC ^ IBS) - (CRC ^ IBD ^ IBS)
CRC_specific <- CRC_case_minus_control[!(CRC_case_minus_control$pathway %in% CRC_IBD_IBS_case_minus_control_overlap$pathway),]
CRC_specific <- CRC_specific[!(CRC_specific$pathway %in% CRC_IBD_case_minus_control_overlap$pathway),]
CRC_specific <- CRC_specific[!(CRC_specific$pathway %in% CRC_IBS_case_minus_control_overlap$pathway),]
dim(CRC_specific) 

## sort 
CRC_specific <- CRC_specific[order(CRC_specific$p_adj),]


## IBD_specific = IBD - (CRC ^ IBD) - (IBD ^ IBS) - (IBD ^ IBD ^ IBS)
IBD_specific <- IBD_case_minus_control[!(IBD_case_minus_control$pathway %in% CRC_IBD_IBS_case_minus_control_overlap$pathway),]
IBD_specific <- IBD_specific[!(IBD_specific$pathway %in% CRC_IBD_case_minus_control_overlap$pathway),]
IBD_specific <- IBD_specific[!(IBD_specific$pathway %in% IBD_IBS_case_minus_control_overlap$pathway),]
dim(IBD_specific)

IBD_specific <- IBD_specific[order(IBD_specific$p_adj),]

## IBS_specific = IBS - (CRC ^ IBS) - (IBS ^ IBS) - (IBS ^ IBS ^ IBS)
IBS_specific <- IBS_case_minus_control[!(IBS_case_minus_control$pathway %in% CRC_IBD_IBS_case_minus_control_overlap$pathway),]
IBS_specific <- IBS_specific[!(IBS_specific$pathway %in% CRC_IBS_case_minus_control_overlap$pathway),]
IBS_specific <- IBS_specific[!(IBS_specific$pathway %in% IBD_IBS_case_minus_control_overlap$pathway),]
dim(IBS_specific)


IBS_specific <- IBS_specific[order(IBS_specific$p_adj),]

################ Combine top pathways from each overlapped set #################

## Extract the relevant columns before combining them
CRC_IBD_IBS_top_path <- CRC_IBD_IBS_case_minus_control_overlap[,c(1,2,9,16,17,20,27,34,35,38,45,52,53)]
CRC_IBD_top_path <- CRC_IBD_case_minus_control_overlap[,c(1,2,9,16,17,20,27,34,35)]
CRC_IBS_top_path <- CRC_IBS_case_minus_control_overlap[,c(1,2,9,16,17,20,27,34,35)]
IBD_IBS_top_path <- IBD_IBS_case_minus_control_overlap[,c(1,2,9,16,17,20,27,34,35)]
CRC_top_path <- CRC_specific[,c(3,1,9,16,17)]
IBD_top_path <- IBD_specific[,c(3,1,9,16,17)]
IBS_top_path <- IBS_specific[,c(3,1,9,16,17)]

## Extract the relevant columns before combining them
colnames(CRC_IBD_IBS_top_path)
colnames(CRC_IBD_IBS_top_path) <- gsub("\\.x","_CRC",colnames(CRC_IBD_IBS_top_path))
colnames(CRC_IBD_IBS_top_path) <- gsub("\\.y","_IBD",colnames(CRC_IBD_IBS_top_path))
colnames(CRC_IBD_IBS_top_path)[10:13] <- paste(colnames(CRC_IBD_IBS_top_path)[10:13],"_IBS",sep ="")
colnames(CRC_IBD_IBS_top_path)

colnames(CRC_IBD_top_path)
colnames(CRC_IBD_top_path) <- gsub("\\.x","_CRC",colnames(CRC_IBD_top_path))
colnames(CRC_IBD_top_path) <- gsub("\\.y","_IBD",colnames(CRC_IBD_top_path))
colnames(CRC_IBD_top_path)

colnames(CRC_IBS_top_path)
colnames(CRC_IBS_top_path) <- gsub("\\.x","_CRC",colnames(CRC_IBS_top_path))
colnames(CRC_IBS_top_path) <- gsub("\\.y","_IBS",colnames(CRC_IBS_top_path))
colnames(CRC_IBS_top_path)

colnames(IBD_IBS_top_path)
colnames(IBD_IBS_top_path) <- gsub("\\.x","_IBD",colnames(IBD_IBS_top_path))
colnames(IBD_IBS_top_path) <- gsub("\\.y","_IBS",colnames(IBD_IBS_top_path))
colnames(IBD_IBS_top_path)


colnames(CRC_top_path) 
colnames(CRC_top_path)[2:5] <- paste(colnames(CRC_top_path)[2:5],"_CRC",sep ="")
colnames(CRC_top_path)

colnames(IBD_top_path) 
colnames(IBD_top_path)[2:5] <- paste(colnames(IBD_top_path)[2:5],"_IBD",sep ="")
colnames(IBD_top_path)

colnames(IBS_top_path) 
colnames(IBS_top_path)[2:5] <- paste(colnames(IBS_top_path)[2:5],"_IBS",sep ="")
colnames(IBS_top_path)

combined_top_path <- rbind(CRC_IBD_IBS_top_path,CRC_IBD_top_path,CRC_IBS_top_path,IBD_IBS_top_path,
                           CRC_top_path,IBD_top_path,IBS_top_path) 

## Add additional columns to dataframes to be able to rbind.
CRC_IBD_top_path[setdiff(names(CRC_IBD_IBS_top_path), names(CRC_IBD_top_path))] <- 0
CRC_IBS_top_path[setdiff(names(CRC_IBD_IBS_top_path), names(CRC_IBS_top_path))] <- 0
IBD_IBS_top_path[setdiff(names(CRC_IBD_IBS_top_path), names(IBD_IBS_top_path))] <- 0
CRC_top_path[setdiff(names(CRC_IBD_IBS_top_path), names(CRC_top_path))] <- 0
IBD_top_path[setdiff(names(CRC_IBD_IBS_top_path), names(IBD_top_path))] <- 0
IBS_top_path[setdiff(names(CRC_IBD_IBS_top_path), names(IBS_top_path))] <- 0

combined_top_path <- rbind(CRC_IBD_IBS_top_path,CRC_IBD_top_path,CRC_IBS_top_path,IBD_IBS_top_path,
                           CRC_top_path,IBD_top_path,IBS_top_path) 

