## Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu

## This script identifies case-specific associations and features for CRC, IBD and IBS
## by factoring out control associations/features from case associations/features.

## Initial setup 
rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
setwd("/path/to/working_directory")

############## Input datasets ##############

## CRC
## case
dataset <- "CRC/case"
filename <- "CRC_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("./data/analysis/", dataset, "/output_lasso_hdi_loocv/combined_output/", filename)
case <- data.frame(fread(filepath, sep="\t", head=T), row.names = 1)
dim(case)

## control 
dataset <- "CRC/control"
filename <- "CRC_control_gene_taxa_FDR_0.1.txt"
filepath <- paste0("./data/analysis/", dataset, "/output_lasso_hdi_loocv/combined_output/", filename)
control <- data.frame(fread(filepath, sep="\t", head=T), row.names = 1)
dim(control)

## IBD
## case
dataset <- "IBD/case"
filename <- "IBD_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("./data/analysis/", dataset, "/output_lasso_hdi_loocv/combined_output/", filename)
case <- data.frame(fread(filepath, sep="\t", head=T), row.names = 1)
dim(case)

## control 
dataset <- "IBD/control"
filename <- "IBD_control_gene_taxa_FDR_0.1.txt"
filepath <- paste0("./data/analysis/", dataset, "/output_lasso_hdi_loocv/combined_output/", filename)
control <- data.frame(fread(filepath, sep="\t", head=T), row.names = 1)
dim(control)

## IBS
## case
dataset <- "IBS/case"
filename <- "IBS_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("./data/analysis/", dataset, "/output_lasso_hdi_loocv/combined_output/", filename)
case <- data.frame(fread(filepath, sep="\t", head=T), row.names = 1)
dim(case)

## control 
dataset <- "IBS/control"
filename <- "IBS_control_gene_taxa_FDR_0.1.txt"
filepath <- paste0("./data/analysis/", dataset, "/output_lasso_hdi_loocv/combined_output/", filename)
control <- data.frame(fread(filepath, sep="\t", head=T), row.names = 1)
dim(control)

############# Compute overlaps (FDR < 0.1) ###########

## overlapping associations
case_control_overlap_assoc <- merge(case, control, by = c("gene", "taxa"))
dim(case_control_overlap_assoc) 

## overlapping genes
case_control_overlap_genes <- intersect(case$gene, control$gene)
length(case_control_overlap_genes) 
case_control_overlap_genes

## overlapping taxa
case_control_overlap_taxa <- intersect(case$taxa, control$taxa)
length(case_control_overlap_taxa)
case_control_overlap_taxa

## write only case-specific genes to a file
case_specific_genes <- setdiff(case$gene, control$gene)
length(case_specific_genes)

dataset <- "disease/case"
filename <- "case_specific_genes_assoc_w_taxa.txt"
filepath <- paste0("./data/analysis/", dataset, "/output_lasso_hdi_loocv/combined_output/", filename)
write(case_specific_genes, file = filepath, sep="\n")



