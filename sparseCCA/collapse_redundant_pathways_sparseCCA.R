## Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu

## This script computes similarity metrics between pathways to identify and collapse redundant pathways

## Initial setup 
rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
setwd("/path/to/working_directory")

library(msigdbr)
library(reshape2)
library(igraph)


####### Functions ##########
## Function to compute overlap coefficient
compute_overlap_coeff <- function(selected_paths, pathway_DB){
  
  
  ## create a matrix using input vector of selected pathways
  overlap_matrix <- matrix(0, nrow = length(selected_paths), ncol = length(selected_paths))
  rownames(overlap_matrix) <- selected_paths
  colnames(overlap_matrix) <- selected_paths
  
  ## iterate over matrix to compute pairwise overlap coeff
  for( i in 1:nrow(overlap_matrix)){
    for(j in 1:ncol(overlap_matrix)){
      ## get pairs of pathways
      path_i <- rownames(overlap_matrix)[i]
      path_j <- colnames(overlap_matrix)[j]
      
      ## genes in these pathway
      path_i_gene_set <- pathway_DB[pathway_DB$gs_name == path_i,]$human_gene_symbol
      length(path_i_gene_set) 
      path_j_gene_set <- pathway_DB[pathway_DB$gs_name == path_j,]$human_gene_symbol
      length(path_j_gene_set) 
      
      ## compute overlap coefficient between these two gene sets
      ## Overlap coefficient is defined as: the size of the intersection divided by the smaller of the size of the two sets.
      ## overlap(X,Y) = |X ^ Y|/min(|X|,|Y|)
      path_i_path_j_intersect <- length(intersect(path_i_gene_set,path_j_gene_set)) 
      overlap_coeff <- path_i_path_j_intersect/min(length(path_i_gene_set),length(path_j_gene_set))
      overlap_matrix[i,j] <- overlap_coeff
      
      
    }
  }
  
  return(overlap_matrix)
}

## Use avg log(FDR) across diseases to pick the one with lowest value (most significant pathway)
pick_representative_path_alt <- function(connected_component,pathway_set){
  
  path_avg_FDR <- c()
  for(k in 1:length(connected_component)){
    path <- pathway_set[pathway_set$pathway == connected_component[k],]
    
  }
  return(connected_component[which.min(path_avg_FDR)])
}

########## Read in pathway sets to collapse ##########

## set dirpath
dataset <- "CRC_IBD_IBS"
overlap_dir <- "/output_sparseCCA/kegg_pid/overlap_DE_FDR_cutoff/"
filter <- "25_300_5"

## CRC, IBD, IBS
filename <- "CRC_IBD_IBS_overlap_FDR_0.1.txt"
filepath <- paste0("./data/analysis/",dataset,overlap_dir,
                   "/",filter,"/",filename)
CRC_IBD_IBS_case_minus_control_overlap <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
dim(CRC_IBD_IBS_case_minus_control_overlap)


## CRC, IBD
filename <- "CRC_IBD_overlap_FDR_0.1.txt"
filepath <- paste0("./data/analysis/",dataset,overlap_dir,
                   "/",filter,"/",filename)
CRC_IBD_case_minus_control_overlap <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
dim(CRC_IBD_case_minus_control_overlap) 


## CRC, IBS
filename <- "CRC_IBS_overlap_FDR_0.1.txt"
filepath <- paste0("./data/analysis/",dataset,overlap_dir,
                   "/",filter,"/",filename)
CRC_IBS_case_minus_control_overlap <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
dim(CRC_IBS_case_minus_control_overlap) 


## IBD, IBS
filename <- "IBD_IBS_overlap_FDR_0.1.txt"
filepath <- paste0("./data/analysis/",dataset,overlap_dir,
                   "/",filter,"/",filename)
IBD_IBS_case_minus_control_overlap <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
dim(IBD_IBS_case_minus_control_overlap) 


## CRC-specific
filename <- "CRC_specific_FDR_0.1.txt"
filepath <- paste0("./data/analysis/",dataset,overlap_dir,
                   "/",filter,"/",filename)
CRC_specific <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
dim(CRC_specific) 


## IBD-specific
filename <- "IBD_specific_FDR_0.1.txt"
filepath <- paste0("./data/analysis/",dataset,overlap_dir,
                   "/",filter,"/",filename)
IBD_specific <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
dim(IBD_specific) 


## IBS-specific
filename <- "IBS_specific_FDR_0.1.txt"
filepath <- paste0("./data/analysis/",dataset,overlap_dir,
                   "/",filter,"/",filename)
IBS_specific <- read.table(filepath,sep="\t",head=T, check.names = F, stringsAsFactors = F)
dim(IBS_specific)


########## Fetch gene sets from msigdb ##############
msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
msigdb_C2 <- as.data.frame(msigdb_C2)
msigdb_C2_CP <- msigdb_C2[msigdb_C2$gs_subcat!= "CGP",]
msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
length(msigdb_C2_CP_unique_pathways)

## Pathway DB of interest: KEGG and PID
path_DB <- c("CP:KEGG","CP:PID")
msigdb_C2_CP <- msigdb_C2_CP[(msigdb_C2_CP$gs_subcat %in% path_DB), ]
table(msigdb_C2_CP$gs_subcat)
# CP:KEGG      CP:PID 
# 12783        8058  

## unique pathways in DBs of interest
msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
length(msigdb_C2_CP_unique_pathways) 
#382

########## Compute overlaps for CRC^IBD^IBS ##############
## (CRC ^ IBD ^ IBS)

## set pathway DB
pathway_DB <- msigdb_C2_CP

selected_paths <- CRC_IBD_IBS_case_minus_control_overlap$pathway
## set no. of top path to keep
n_top_path <- 15
if(length(selected_paths) > n_top_path) {
  selected_paths <- selected_paths[1:n_top_path]
}
overlap_CRC_IBD_IBS <- compute_overlap_coeff(selected_paths,pathway_DB)

## Convert matrix to long form so as to only keep upper traingle of matrix
## This is to avoid duplicate elements in the melted output
overlap_CRC_IBD_IBS_long <- melt(replace(overlap_CRC_IBD_IBS, lower.tri(overlap_CRC_IBD_IBS, TRUE), NA), na.rm = TRUE)
dim(overlap_CRC_IBD_IBS_long)

## Identify pairs of pathways with similarity > cut-off
similarity_cutoff <- 0.5
overlap_CRC_IBD_IBS_long_0.3 <- overlap_CRC_IBD_IBS_long[(overlap_CRC_IBD_IBS_long$value > similarity_cutoff),]
dim(overlap_CRC_IBD_IBS_long_0.3)

## Identify connected components from the filtered pairs.
## 1. prepare relevant df
data <- overlap_CRC_IBD_IBS_long_0.3[,1:2]

## 2. Build graph from df
g <- graph.data.frame(d = overlap_CRC_IBD_IBS_long_0.3[,1:2], directed = FALSE)

## 3. partition the graph into disjoint sets
connected_components <- split(V(g)$name, clusters(g)$membership)
connected_components

## Pick a representative pathway
## Compute average log10(FDR) values across the diseases for each pathway within a connected component. 
pathway_set <- CRC_IBD_IBS_case_minus_control_overlap
representative_path <- list()
for(i in 1:length(connected_components)){
  
  representative_path[[i]] <- pick_representative_path_alt(connected_components[[i]],pathway_set)
}
representative_path

## Repeat this process for CRC ^ IBD, CRC ^ IBS, IBD ^ IBS, CRC-specific, IBD-specific, and IBS-specific pathways