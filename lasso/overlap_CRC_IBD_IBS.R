## Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu

## This script checks for overlapping genes, taxa and gene-taxa associations 
## across lasso runs on CRC, IBS and IBD

## Initial setup 
rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
library(UpSetR)
library(VennDiagram)
library(RColorBrewer)
library(msigdbr)
library(data.table)
library(reshape2)
library(ggplot2)

setwd("/path/to/working_directory")

################ Function ##############
msigdb_enrichment <- function(pathway_DB, pathways, background_genes, genes_of_interest){
  
  ## set min and max genes to allow pathways for testing
  min_genes <- 25
  max_genes <- 300
  overlap_genes <- 3
  
  enrichment_list <- list()
  for(i in 1:length(pathways)){
    
    
    pathway <- pathways[i]
    
    ## genes in this pathway
    pathway_gene_set <- pathway_DB[pathway_DB$gs_name == pathway,]$human_gene_symbol
    length(pathway_gene_set) 
    ## If the criteria for min and max #genes in a given pathway is not satified, 
    ## skip testing the current pathway
    if(length(pathway_gene_set) < min_genes || length(pathway_gene_set) > max_genes) next
    
    ## The contingency table
    ## Inspired by: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
    ##                    genes_of_interest genes_NOT_of_interest  Total
    ## In_pathway               x                 m - x              m
    ## Not_in_pathway         k - x             n - (k - x)          n
    ## Total                    k                 m + n - k          m + n
    
    ## m, #overlapping genes in this pathway and background genes for IBS
    m <- length(intersect(background_genes,pathway_gene_set))
    m ## 11
    
    ## n, #genes in background but not in pathway 
    n <- length(setdiff(background_genes,pathway_gene_set))
    n #12502
    
    ## x, #genes of interest in pathway 
    x <- length(intersect(pathway_gene_set,genes_of_interest))
    x ## 1
    ## If the overlap between genes of interest and the pathway is less than the overlap cut-off, 
    ## then skip testing the current pathway
    if(x < overlap_genes) next 
    
    # if(x > 0){
    #   gene_names_in_pathway = noquote(paste(intersect(pathway_gene_set,genes_of_interest), collapse = ","))
    # } else {
    #   gene_names_in_pathway = ""
    # }
    ## Extract list of genes in the genes of interest that are included in the pathway. 
    gene_names_in_pathway = noquote(paste(intersect(pathway_gene_set,genes_of_interest), collapse = ","))
    
    ## k, total #genes in genes list of interest
    k <- length(genes_of_interest)
    k #839
    
    ## Build the contigency table
    contingency_table <- matrix(c(x, k-x, m-x, n-(k-x)), #matrix is filled along the column by default
                                nrow = 2, ncol = 2, 
                                dimnames = list(c("In_pathway","Not_in_pathway"),
                                                c("Genes_of_interest","Genes_NOT_of_interest"))
    )
    
    contingency_table
    
    fisher_result <- fisher.test( contingency_table, alternative = "greater")
    
    ## save details in a dataframe
    enrichment_result_df <- data.frame( pathway = pathway,
                                        BG_genes = length(background_genes),
                                        genes_in_pathway = length(pathway_gene_set),
                                        genes_in_path_and_BG = m,
                                        genes_of_interest = k,
                                        genes_of_interest_in_pathway = x,
                                        gene_names_in_pathway = gene_names_in_pathway,
                                        ## fill in contingency table entries: x, k-x, m-x, n-(k-x) for z-score computation.
                                        cont_n1 = x,
                                        cont_n2 = k-x,
                                        cont_n3 = m-x,
                                        cont_n4 = n-(k-x),
                                        CI_95 = paste0(signif(fisher_result$conf.int[1],5),"-",signif(fisher_result$conf.int[2],5)),
                                        odds_ratio = unname(fisher_result$estimate),
                                        p_val = fisher_result$p.value
    )
    
    enrichment_list[[i]] <- enrichment_result_df
    
  }
  
  return(enrichment_list)
  
}
pick_taxaname <- function(microbe_name){
  
  # microbe_name <- gsub("(NA\\>)+","",microbe_name) # \\> to match end of word, + for once or more.
  # microbe_name <- gsub("(;)+$","",microbe_name,perl=T) # \\> to match end of word, + for once or more.
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
}
compute_spearman_corr <- function(gene_taxa_assoc, genes, microbes){
  
  select.genes <- as.character(unique(gene_taxa_assoc$gene)) 
  select.taxa <- as.character(unique(gene_taxa_assoc$taxa)) 
  
  select.gene.expr <- genes[,which(colnames(genes) %in% select.genes)]
  select.taxa.abund <- microbes[,which(colnames(microbes) %in% select.taxa)]
  
  spear_rho <- cor(select.gene.expr, select.taxa.abund, method = "spearman")
  
  return(spear_rho)
}

################## Read in case-specific lasso output for each disease #############

dataset <- "CRC/case"
filename <- "CRC_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("data/analysis/",dataset,"/output_lasso_hdi_loocv/combined_output/",filename)
CRC.lasso.stabsel <- read.table(filepath,sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
dim(CRC.lasso.stabsel) 

dataset <- "IBD/case"
filename <- "IBD_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("data/analysis/",dataset,"/output_lasso_hdi_loocv/combined_output/",filename)
IBD.lasso.stabsel <- read.table(filepath,sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
dim(IBD.lasso.stabsel)

dataset <- "IBS/case"
filename <- "IBS_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("data/analysis/",dataset,"/output_lasso_hdi_loocv/combined_output/",filename)
IBS.lasso.stabsel <- read.table(filepath,sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
dim(IBS.lasso.stabsel)

################ Input gene expr and taxa abundance data for each disease (case) ##########

## CRC
dataset <- "CRC"
case_control <- "case"
## input genes 
genes.file <-"CRC_genes.txt"
input.genes.table <- paste0("./data/clean/",dataset,"/gene_expr/",genes.file) # samples x genes
CRC_case_genes <- data.frame(fread(input.genes.table,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F) ## check.names = F doesn't change gene leabels like "ERVFRD-1" to "ERVFRD.1"
dim(CRC_case_genes)

## input microbiome
taxa.file <- "CRC_taxa.txt"
input.microbiome.table <- paste0("./data/clean/",dataset,"/taxa_abnd/",taxa.file) # microbiome x samples
CRC_case_microbes <- read.table(input.microbiome.table,sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
dim(CRC_case_microbes)

## IBD
dataset <- "IBD"
case_control <- "case"
## input genes 
genes.file <-"IBD_genes.txt"
input.genes.table <- paste0("./data/clean/",dataset,"/gene_expr/",genes.file) # samples x genes
IBD_case_genes <- data.frame(fread(input.genes.table,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
dim(IBD_case_genes)

## input microbiome
taxa.file <- "IBD_taxa.txt"
input.microbiome.table <- paste0("./data/clean/",dataset,"/taxa_abnd/",taxa.file) # microbiome x samples
IBD_case_microbes <- read.table(input.microbiome.table,sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
dim(IBD_case_microbes) 

## IBS
dataset <- "IBS"
case_control <- "case"
## input genes 
genes.file <-"IBS_genes.txt"
input.genes.table <- paste0("./data/clean/",dataset,"/gene_expr/",genes.file) # samples x genes
IBS_case_genes <- data.frame(fread(input.genes.table,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
dim(IBS_case_genes)

## input microbiome
taxa.file <- "IBS_taxa.txt"
input.microbiome.table <- paste0("./data/clean/",dataset,"/taxa_abnd/",taxa.file) # microbiome x samples
IBS_case_microbes <- read.table(input.microbiome.table,sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
dim(IBS_case_microbes)

################ Compute spearman correlation for each datasets gene-taxa pairs ###########
## This will be used for visualizing networks

# CRC
## set assocation to plot
gene_taxa_assoc <- CRC.lasso.stabsel
dim(gene_taxa_assoc) 
genes <- CRC_case_genes 
microbes <- CRC_case_microbes

## Compute spearman correlation 
spear_rho <- compute_spearman_corr(gene_taxa_assoc, genes, microbes)
dim(spear_rho)
any(is.na(spear_rho)) # F


## melt spearman output to long format
spear_rho_melt <- melt(spear_rho)
dim(spear_rho_melt)
colnames(spear_rho_melt) <- c("gene","taxa", "spear_rho")

## subset the correlation to gene-taxa pairs in the lasso output
CRC.lasso.stabsel.spear <- merge(CRC.lasso.stabsel,spear_rho_melt, by = c("gene","taxa"))
dim(CRC.lasso.stabsel.spear)
## sort by BH_global
CRC.lasso.stabsel.spear <- CRC.lasso.stabsel.spear[order(CRC.lasso.stabsel.spear$BH_global),]
CRC.lasso.stabsel <- CRC.lasso.stabsel.spear


## IBD
## set assocation to plot
gene_taxa_assoc <- IBD.lasso.stabsel
dim(gene_taxa_assoc) 
genes <- IBD_case_genes 
microbes <- IBD_case_microbes

## Compute spearman correlation 
spear_rho <- compute_spearman_corr(gene_taxa_assoc, genes, microbes)
dim(spear_rho)
any(is.na(spear_rho)) # F


## melt spearman output to long format
spear_rho_melt <- melt(spear_rho)
dim(spear_rho_melt)
colnames(spear_rho_melt) <- c("gene","taxa", "spear_rho")

## subset the correlation to gene-taxa pairs in the lasso+stabsel output
IBD.lasso.stabsel.spear <- merge(IBD.lasso.stabsel,spear_rho_melt, by = c("gene","taxa"))
dim(IBD.lasso.stabsel.spear)
## sort by BH_global
IBD.lasso.stabsel.spear <- IBD.lasso.stabsel.spear[order(IBD.lasso.stabsel.spear$BH_global),]
IBD.lasso.stabsel <- IBD.lasso.stabsel.spear

## IBS
## set assocation to plot
gene_taxa_assoc <- IBS.lasso.stabsel
dim(gene_taxa_assoc) 
genes <- IBS_case_genes 
microbes <- IBS_case_microbes

## Compute spearman correlation 
spear_rho <- compute_spearman_corr(gene_taxa_assoc, genes, microbes)
dim(spear_rho) 
any(is.na(spear_rho))

## melt spearman output to long format
spear_rho_melt <- melt(spear_rho)
dim(spear_rho_melt)
colnames(spear_rho_melt) <- c("gene","taxa", "spear_rho")

## subset the correlation to gene-taxa pairs in the lasso+stabsel output
IBS.lasso.stabsel.spear <- merge(IBS.lasso.stabsel,spear_rho_melt, by = c("gene","taxa"))
dim(IBS.lasso.stabsel.spear)
## sort by BH_global
IBS.lasso.stabsel.spear <- IBS.lasso.stabsel.spear[order(IBS.lasso.stabsel.spear$BH_global),]
IBS.lasso.stabsel <- IBS.lasso.stabsel.spear

################ Test overlaps pair-wise and across three diseases ################
## Unique genes per dataset
CRC.genes <- unique(CRC.lasso.stabsel$gene); length(CRC.genes) 

IBD.genes <- unique(IBD.lasso.stabsel$gene); length(IBD.genes) 

IBS.genes <- unique(IBS.lasso.stabsel$gene); length(IBS.genes) 


## overlapping gene sets
## CRC ^ IBD ^ IBS
CRC.IBD.IBS.genes <- Reduce(intersect, list(CRC.genes,IBD.genes,IBS.genes))
length(CRC.IBD.IBS.genes)
CRC.IBD.IBS.genes

## CRC ^ IBD = CRC ^ IBD - CRC ^ IBD ^ IBS
CRC.IBD.genes.overlap <- setdiff(intersect(CRC.genes,IBD.genes), CRC.IBD.IBS.genes)
length(CRC.IBD.genes.overlap)

## CRC ^ IBS
CRC.IBS.genes.overlap <- setdiff(intersect(CRC.genes, IBS.genes), CRC.IBD.IBS.genes)
length(CRC.IBS.genes.overlap) 
IBD.IBS.genes.overlap <- setdiff(intersect(IBD.genes, IBS.genes), CRC.IBD.IBS.genes)
length(IBD.IBS.genes.overlap) 


## Unique taxa per dataset
CRC.taxa <- unique(CRC.lasso.stabsel$taxa); length(CRC.taxa)

IBD.taxa <- unique(IBD.lasso.stabsel$taxa); length(IBD.taxa) 

IBS.taxa <- unique(IBS.lasso.stabsel$taxa); length(IBS.taxa) 

## overlapping taxa sets
## CRC ^ IBD ^ IBS
CRC.IBD.IBS.taxa <- Reduce(intersect, list(CRC.taxa,IBD.taxa,IBS.taxa))
CRC.IBD.IBS.taxa

## CRC ^ IBD = CRC ^ IBD - CRC ^ iBD ^ IBS
CRC.IBD.taxa.overlap <- setdiff(intersect(CRC.taxa,IBD.taxa), CRC.IBD.IBS.taxa)
length(CRC.IBD.taxa.overlap) 

## CRC ^ IBS
CRC.IBS.taxa.overlap <-  setdiff(intersect(CRC.taxa, IBS.taxa), CRC.IBD.IBS.taxa)
length(CRC.IBS.taxa.overlap)

IBD.IBS.taxa.overlap <- setdiff(intersect(IBD.taxa, IBS.taxa), CRC.IBD.IBS.taxa)
length(IBD.IBS.taxa.overlap)

## Overlapping gene-taxa interactions 
## CRC ^ IBD ^ IBS
CRC.IBD.IBS.assoc.overlap <- Reduce(function(x,y) merge(x=x, y=y, by = c("gene","taxa")), list(CRC.lasso.stabsel,IBD.lasso.stabsel, IBS.lasso.stabsel))
dim(CRC.IBD.IBS.assoc.overlap)

## CRC ^ IBD
CRC.IBD.assoc.overlap <- merge(CRC.lasso.stabsel,IBD.lasso.stabsel, by = c("gene","taxa"))
dim(CRC.IBD.assoc.overlap) 
CRC.IBD.assoc.overlap

##IBD ^ IBS
IBD.IBS.assoc.overlap <- merge(IBD.lasso.stabsel,IBS.lasso.stabsel, by = c("gene","taxa"))
dim(IBD.IBS.assoc.overlap)

## CRC ^ IBS
CRC.IBS.assoc.overlap <- merge(CRC.lasso.stabsel,IBS.lasso.stabsel, by = c("gene","taxa"))
dim(CRC.IBS.assoc.overlap)
################# Perform enrichment for overlapping gene sets using msigdb ##############

## Create background gene set --  union of CRC, IBD, and IBS case genes will be background. 
## CRC case genes
filename <- "./data/clean/CRC/gene_expr/CRC_genes.txt"
genes <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
genes <- t(genes)
dim(genes)
CRC.case.background.genes  <- rownames(genes)
CRC.case.background.genes <- CRC.case.background.genes[order(CRC.case.background.genes)]
length(CRC.case.background.genes)

## IBD case genes
filename <- "./data/clean/IBD/gene_expr/IBD_genes.txt"
genes <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
genes <- t(genes)
dim(genes) 
IBD.case.background.genes  <- rownames(genes)
IBD.case.background.genes <- IBD.case.background.genes[order(IBD.case.background.genes)]
length(IBD.case.background.genes) 

## IBS case genes
filename <- "./data/clean/IBS/gene_expr/IBS_genes.txt"
genes <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
genes <- t(genes)
dim(genes) 
IBS.case.background.genes  <- rownames(genes)
IBS.case.background.genes <- IBS.case.background.genes[order(IBS.case.background.genes)]
length(IBS.case.background.genes) 

## Background genes  = union of CRC, IBD and IBS genes
background_genes <- Reduce(union, list(CRC.case.background.genes,IBD.case.background.genes,IBS.case.background.genes))
length(background_genes)

## Get MSIGDB pathways of interest
msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
class(msigdb_C2) #[1] "tbl_df"     "tbl"        "data.frame"
msigdb_C2 <- as.data.frame(msigdb_C2)
table(msigdb_C2$gs_subcat)

msigdb_C2_CP <- msigdb_C2[msigdb_C2$gs_subcat!= "CGP",]
table(msigdb_C2_CP$gs_subcat)

dim(msigdb_C2_CP)
length(unique(msigdb_C2_CP$gs_name))

## Only keep pathway DBs that we want to test gene sets against 
table(msigdb_C2_CP$gs_subcat)
path_DB <- c("CP:KEGG","CP:PID")
msigdb_C2_CP <- msigdb_C2_CP[(msigdb_C2_CP$gs_subcat %in% path_DB), ]
table(msigdb_C2_CP$gs_subcat)

## unique pathways in DBs of interest
msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
length(msigdb_C2_CP_unique_pathways) 

## Perform enrichment
pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways
background_genes <- background_genes


## Check different overlp gene sets
## CRC ^ IBD includes some significant pathways
genes_of_interest <- CRC.IBD.genes.overlap

##Call enrichment function
overlap_geneset_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
## convert the list to dataframe
msigdb_pathways_df <- do.call(rbind.data.frame, overlap_geneset_pathways)
dim(msigdb_pathways_df)
##sort by pval
msigdb_pathways_df <- msigdb_pathways_df[order(msigdb_pathways_df$p_val),]
msigdb_pathways_df$p_adj <- p.adjust(msigdb_pathways_df$p_val, method = "BH")
##how many pathways at FDR < 0.1
length(which(msigdb_pathways_df$p_adj < 0.1))

################# Venn Diagram to visualize overlapping genes and taxa ############
listGenes <- list(CRC = unique(as.character(CRC.genes)), IBD = unique(as.character(IBD.genes)), 
                  IBS = unique((as.character(IBS.genes))))

listTaxa <- list(CRC = as.character(CRC.taxa), IBD = as.character(IBD.taxa), 
                 IBS = as.character(IBS.taxa))

dataset <- "CRC_IBD_IBS"


## V3
crc.col.fill <- "#dfb5ff" 
ibd.col.fill <- "#a7f8fa"
ibs.col.fill <- "#ebed79"
myCol <- c(crc.col.fill,ibd.col.fill,ibs.col.fill)

dataset <- "CRC_IBD_IBS"
filename <- "CRC_IBD_IBS_overlapping_genes_venn6.pdf" ## svg version needed for inkscape (vector graphics)

genes_venn<- venn.diagram(
  x = listGenes,
  category.names = c("CRC" , "IBD" , "IBS"),
  # filename = paste0("./data/analysis/",dataset,"/plots/lasso_overlap/",filename),
  filename = NULL, ## for pdf output
  output=TRUE, 
  
  # Output features
  # imagetype="svg" ,
  height = 480 , 
  width = 480 , 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 2,
  # fontface = "bold",
  # fontfamily = "serif",
  
  # Set names
  cat.cex = 2,
  # cat.fontface = "bold",
  cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 180),
  cat.pos = c(-27, 27, 180),
  cat.dist = c(0.055, 0.055, 0.055),
  # cat.fontfamily = "serif",
  rotation = 1
  
)
ggsave(genes_venn, file=paste0("./data/analysis/",dataset,"/plots/lasso_overlap/",filename), device = "pdf")

## Repeat for taxa
dataset <- "CRC_IBD_IBS"
filename <- "CRC_IBD_IBS_overlapping_taxa_venn6.pdf"

taxa_venn <- venn.diagram(
  x = listTaxa,
  category.names = c("CRC" , "IBD" , "IBS"),
  # filename = paste0("./data/analysis/",dataset,"/plots/lasso_overlap/",filename),
  filename = NULL, ## for pdf output
  output=TRUE, 
  
  # Output features
  # imagetype="svg" , ## for pdf output
  height = 480 , 
  width = 480 , 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 2,
  # fontface = "bold",
  # fontfamily = "sans",
  
  # Set names
  cat.cex = 2,
  # cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 180),
  cat.dist = c(0.055, 0.055, 0.055),
  # cat.fontfamily = "sans",
  rotation = 1
)
ggsave(taxa_venn, file=paste0("./data/analysis/",dataset,"/plots/lasso_overlap/",filename), device = "pdf")


################# Identify disease-specific genes for overlapping taxa ################

## CRC, IBD
CRC.IBD.taxa.overlap

CRC.genes.in.CRC.IBD.taxa.overlap <- CRC.lasso.stabsel[(CRC.lasso.stabsel$taxa %in% CRC.IBD.taxa.overlap),]$gene
length(CRC.genes.in.CRC.IBD.taxa.overlap) 
IBD.genes.in.CRC.IBD.taxa.overlap <- IBD.lasso.stabsel[(IBD.lasso.stabsel$taxa %in% CRC.IBD.taxa.overlap),]$gene
length(IBD.genes.in.CRC.IBD.taxa.overlap) 
intersect(CRC.genes.in.CRC.IBD.taxa.overlap,IBD.genes.in.CRC.IBD.taxa.overlap) 
 

## CRC, IBS
CRC.IBS.taxa.overlap

CRC.genes.in.CRC.IBS.taxa.overlap <- CRC.lasso.stabsel[(CRC.lasso.stabsel$taxa %in% CRC.IBS.taxa.overlap),]$gene
length(CRC.genes.in.CRC.IBS.taxa.overlap) 
CRC.genes.in.CRC.IBS.taxa.overlap
IBS.genes.in.CRC.IBS.taxa.overlap <- IBS.lasso.stabsel[(IBS.lasso.stabsel$taxa %in% CRC.IBS.taxa.overlap),]$gene
length(IBS.genes.in.CRC.IBS.taxa.overlap) 
IBS.genes.in.CRC.IBS.taxa.overlap
## any overlap between these genes 
intersect(CRC.genes.in.CRC.IBS.taxa.overlap,IBS.genes.in.CRC.IBS.taxa.overlap) 

## IBD, IBS
IBD.IBS.taxa.overlap

IBD.genes.in.IBD.IBS.taxa.overlap <- IBD.lasso.stabsel[(IBD.lasso.stabsel$taxa %in% IBD.IBS.taxa.overlap),]$gene
length(IBD.genes.in.IBD.IBS.taxa.overlap) 
IBS.genes.in.IBD.IBS.taxa.overlap <- IBS.lasso.stabsel[(IBS.lasso.stabsel$taxa %in% IBD.IBS.taxa.overlap),]$gene
length(IBS.genes.in.IBD.IBS.taxa.overlap) 
## any overlap between these genes 
intersect(IBD.genes.in.IBD.IBS.taxa.overlap,IBS.genes.in.IBD.IBS.taxa.overlap) 

## CRC ^ IBD ^ IBS taxa overlap
CRC.IBD.IBS.taxa

CRC.genes.in.CRC.IBD.IBS.taxa.overlap <- CRC.lasso.stabsel[(CRC.lasso.stabsel$taxa %in% CRC.IBD.IBS.taxa),]$gene
length(CRC.genes.in.CRC.IBD.IBS.taxa.overlap) 
IBD.genes.in.CRC.IBD.IBS.taxa.overlap <- IBD.lasso.stabsel[(IBD.lasso.stabsel$taxa %in% CRC.IBD.IBS.taxa),]$gene
length(IBD.genes.in.CRC.IBD.IBS.taxa.overlap) 
IBS.genes.in.CRC.IBD.IBS.taxa.overlap <- IBS.lasso.stabsel[(IBS.lasso.stabsel$taxa %in% CRC.IBD.IBS.taxa),]$gene
length(IBS.genes.in.CRC.IBD.IBS.taxa.overlap)

################# Identify interactions for (CRC ^ IBD ^ IBS) taxa and plot them #############

## CRC ^ IBD ^ IBS
length(CRC.IBD.IBS.taxa)
## Get CRC assoc for these taxa
CRC.gene.taxa.for.CRC.IBS.IBD <- CRC.lasso.stabsel[(CRC.lasso.stabsel$gene %in% CRC.genes.in.CRC.IBD.IBS.taxa.overlap & 
                                                      CRC.lasso.stabsel$taxa %in% CRC.IBD.IBS.taxa),]

CRC.gene.taxa.for.CRC.IBS.IBD$dataset <- "CRC"
dim(CRC.gene.taxa.for.CRC.IBS.IBD)

## if number of assoc is > 10, pick top 10 assoc., else keep all assoc
gene_taxa_assoc <- CRC.gene.taxa.for.CRC.IBS.IBD
gene_taxa_assoc.filt <- data.frame()
if(dim(gene_taxa_assoc)[1] > 10){
  ## get top interactions which includes all the above shared taxa
  
  while(dim(gene_taxa_assoc.filt)[1] < 10 & dim(gene_taxa_assoc)[1] > 0){
    gene_taxa_assoc.dt <- setDT(gene_taxa_assoc)[order(BH_global), .SD[1], by=taxa]
    gene_taxa_assoc.filt <- rbind(gene_taxa_assoc.filt,as.data.frame(gene_taxa_assoc.dt))
    # dim(gene_taxa_assoc.filt)
    gene_taxa_assoc <- gene_taxa_assoc[ !(gene_taxa_assoc$gene %in% gene_taxa_assoc.filt$gene & 
                                            gene_taxa_assoc$taxa %in% gene_taxa_assoc.filt$taxa),]
    
    
  }
  dim(gene_taxa_assoc.filt)
  ## if greater than 10 rows, subset to pick top 10
  CRC.gene.taxa.for.CRC.IBS.IBD <- gene_taxa_assoc.filt
}

dim(CRC.gene.taxa.for.CRC.IBS.IBD)


## Get IBD assoc for these taxa
IBD.gene.taxa.for.CRC.IBS.IBD <- IBD.lasso.stabsel[(IBD.lasso.stabsel$gene %in% IBD.genes.in.CRC.IBD.IBS.taxa.overlap & 
                                                      IBD.lasso.stabsel$taxa %in% CRC.IBD.IBS.taxa),]

IBD.gene.taxa.for.CRC.IBS.IBD$dataset <- "IBD"
dim(IBD.gene.taxa.for.CRC.IBS.IBD)



## if number of assoc is > 10, pick top 10 assoc., else keep all assoc
gene_taxa_assoc <- IBD.gene.taxa.for.CRC.IBS.IBD
gene_taxa_assoc.filt <- data.frame()
if(dim(gene_taxa_assoc)[1] > 10){
## get top interactions which includes all the above shared taxa
 
  while(dim(gene_taxa_assoc.filt)[1] < 10 & dim(gene_taxa_assoc)[1] > 0){
  gene_taxa_assoc.dt <- setDT(gene_taxa_assoc)[order(BH_global), .SD[1], by=taxa]
  gene_taxa_assoc.filt <- rbind(gene_taxa_assoc.filt,as.data.frame(gene_taxa_assoc.dt))
  # dim(gene_taxa_assoc.filt)
  gene_taxa_assoc <- gene_taxa_assoc[ !(gene_taxa_assoc$gene %in% gene_taxa_assoc.filt$gene & 
                                          gene_taxa_assoc$taxa %in% gene_taxa_assoc.filt$taxa),]
  
  
  }
  dim(gene_taxa_assoc.filt) 
  ## if greater than 10 rows, subset to pick top 10
  IBD.gene.taxa.for.CRC.IBS.IBD <- gene_taxa_assoc.filt
  
}

dim(IBD.gene.taxa.for.CRC.IBS.IBD)

## Get IBS assoc for these taxa
IBS.gene.taxa.for.CRC.IBS.IBD <- IBS.lasso.stabsel[(IBS.lasso.stabsel$gene %in% IBS.genes.in.CRC.IBD.IBS.taxa.overlap & 
                                                      IBS.lasso.stabsel$taxa %in% CRC.IBD.IBS.taxa),]
dim(IBS.gene.taxa.for.CRC.IBS.IBD) 
IBS.gene.taxa.for.CRC.IBS.IBD$dataset <- "IBS"

## get top interactions which includes all the above shared taxa
## if number of assoc is > 10, pick top 10 assoc., else keep all assoc
gene_taxa_assoc <- IBS.gene.taxa.for.CRC.IBS.IBD
gene_taxa_assoc.filt <- data.frame()
if(dim(gene_taxa_assoc)[1] > 10){
  ## get top interactions which includes all the above shared taxa
  
  while(dim(gene_taxa_assoc.filt)[1] < 10 & dim(gene_taxa_assoc)[1] > 0){
    gene_taxa_assoc.dt <- setDT(gene_taxa_assoc)[order(BH_global), .SD[1], by=taxa]
    gene_taxa_assoc.filt <- rbind(gene_taxa_assoc.filt,as.data.frame(gene_taxa_assoc.dt))
    dim(gene_taxa_assoc.filt)
    gene_taxa_assoc <- gene_taxa_assoc[ !(gene_taxa_assoc$gene %in% gene_taxa_assoc.filt$gene & 
                                            gene_taxa_assoc$taxa %in% gene_taxa_assoc.filt$taxa),]
    dim(gene_taxa_assoc) 
    
  }
  dim(gene_taxa_assoc.filt) 
  ## if greater than 10 rows, subset to pick top 10
  IBS.gene.taxa.for.CRC.IBS.IBD <- gene_taxa_assoc.filt
  
}

dim(IBS.gene.taxa.for.CRC.IBS.IBD) 

## combine these assoc. for CRC ^ IBD ^ IBS taxa for Viz.
combined_interactions <- rbind(CRC.gene.taxa.for.CRC.IBS.IBD,
                                  IBD.gene.taxa.for.CRC.IBS.IBD,
                                  IBS.gene.taxa.for.CRC.IBS.IBD
)
dim(combined_interactions) 

combined_interactions$taxa_name <- sapply(as.character(combined_interactions$taxa), pick_taxaname)
combined_interactions$shape_gene <- rep("gene",length(combined_interactions$gene))
combined_interactions$shape_taxa <- rep("taxa",length(combined_interactions$taxa))


combined_interactions$thickness <- abs(combined_interactions$spear_rho)
combined_interactions$sign <- sign(combined_interactions$ci.upper) 

## subset to relevant columns for viz
combined_interactions <- combined_interactions[c(1,2,13:20)]

################# Identify interactions for (CRC ^ IBD) taxa and plot them #############

# CRC ^ IBD
## Extract CRC associations
length(CRC.IBD.taxa.overlap)
CRC.gene.taxa.for.CRC.IBD <- CRC.lasso.stabsel[(CRC.lasso.stabsel$gene %in% CRC.genes.in.CRC.IBD.taxa.overlap & 
                                                  CRC.lasso.stabsel$taxa %in% CRC.IBD.taxa.overlap),]
CRC.gene.taxa.for.CRC.IBD$dataset <- "CRC"
dim(CRC.gene.taxa.for.CRC.IBD) 

### *** This is a special case when we want to add an assoc in CRC which also occurs in IBD ***
RIPK3_Blautia_assoc <- CRC.gene.taxa.for.CRC.IBD[CRC.gene.taxa.for.CRC.IBD$gene == "RIPK3",]


## get top interactions which includes all the above shared taxa
## if number of assoc is > 10, pick top 10 assoc., else keep all assoc
gene_taxa_assoc <- CRC.gene.taxa.for.CRC.IBD
gene_taxa_assoc.filt <- data.frame()
if(dim(gene_taxa_assoc)[1] > 10){
  ## get top interactions which includes all the above shared taxa
  
  while(dim(gene_taxa_assoc.filt)[1] < 10 & dim(gene_taxa_assoc)[1] > 0){
    gene_taxa_assoc.dt <- setDT(gene_taxa_assoc)[order(BH_global), .SD[1], by=taxa]
    gene_taxa_assoc.filt <- rbind(gene_taxa_assoc.filt,as.data.frame(gene_taxa_assoc.dt))
    dim(gene_taxa_assoc.filt)
    gene_taxa_assoc <- gene_taxa_assoc[ !(gene_taxa_assoc$gene %in% gene_taxa_assoc.filt$gene & 
                                            gene_taxa_assoc$taxa %in% gene_taxa_assoc.filt$taxa),]
    dim(gene_taxa_assoc) 
    
  }
  dim(gene_taxa_assoc.filt) 
  ## if greater than 10 rows, subset to pick top 10
  CRC.gene.taxa.for.CRC.IBD <- gene_taxa_assoc.filt
  
}

dim(CRC.gene.taxa.for.CRC.IBD) 
## greater than 10 association, subset to top 10
CRC.gene.taxa.for.CRC.IBD <- CRC.gene.taxa.for.CRC.IBD[1:10,]
dim(CRC.gene.taxa.for.CRC.IBD) 

## Add RIPK3 - Blautia association here
CRC.gene.taxa.for.CRC.IBD <- rbind(CRC.gene.taxa.for.CRC.IBD,RIPK3_Blautia_assoc)
dim(CRC.gene.taxa.for.CRC.IBD) 

## extract IBD assoc.
IBD.gene.taxa.for.CRC.IBD <- IBD.lasso.stabsel[(IBD.lasso.stabsel$gene %in% IBD.genes.in.CRC.IBD.taxa.overlap & 
                                                  IBD.lasso.stabsel$taxa %in% CRC.IBD.taxa.overlap),]
IBD.gene.taxa.for.CRC.IBD$dataset <- "IBD"
dim(IBD.gene.taxa.for.CRC.IBD)

## get top interactions which includes all the above shared taxa
## if number of assoc is > 10, pick top 10 assoc., else keep all assoc
gene_taxa_assoc <- IBD.gene.taxa.for.CRC.IBD
gene_taxa_assoc.filt <- data.frame()
if(dim(gene_taxa_assoc)[1] > 10){
  ## get top interactions which includes all the above shared taxa
   while(dim(gene_taxa_assoc.filt)[1] < 10 & dim(gene_taxa_assoc)[1] > 0){
    gene_taxa_assoc.dt <- setDT(gene_taxa_assoc)[order(BH_global), .SD[1], by=taxa]
    gene_taxa_assoc.filt <- rbind(gene_taxa_assoc.filt,as.data.frame(gene_taxa_assoc.dt))
    # dim(gene_taxa_assoc.filt)
    gene_taxa_assoc <- gene_taxa_assoc[ !(gene_taxa_assoc$gene %in% gene_taxa_assoc.filt$gene & 
                                            gene_taxa_assoc$taxa %in% gene_taxa_assoc.filt$taxa),]
    
    
   }
  
  ## if greater than 10 rows, subset to pick top 10
  IBD.gene.taxa.for.CRC.IBD <- gene_taxa_assoc.filt
  
}
dim(IBD.gene.taxa.for.CRC.IBD) 

## combine these assoc. for CRC ^ IBD 
combined_interactions <- rbind(CRC.gene.taxa.for.CRC.IBD,
                               IBD.gene.taxa.for.CRC.IBD
                               )
dim(combined_interactions) 

combined_interactions$taxa_name <- sapply(as.character(combined_interactions$taxa), pick_taxaname)
combined_interactions$shape_gene <- rep("gene",length(combined_interactions$gene))
combined_interactions$shape_taxa <- rep("taxa",length(combined_interactions$taxa))


combined_interactions$thickness <- abs(combined_interactions$spear_rho)
combined_interactions$sign <- sign(combined_interactions$ci.upper) 

## subset to relevant columns for viz
combined_interactions <- combined_interactions[c(1,2,13:20)]

################# Identify interactions for (CRC ^ IBS) taxa and plot them  #############

## CRC ^ IBS
length(CRC.IBS.taxa.overlap) 
CRC.gene.taxa.for.CRC.IBS <- CRC.lasso.stabsel[(CRC.lasso.stabsel$gene %in% CRC.genes.in.CRC.IBS.taxa.overlap & 
                                                  CRC.lasso.stabsel$taxa %in% CRC.IBS.taxa.overlap),]
dim(CRC.gene.taxa.for.CRC.IBS) 
CRC.gene.taxa.for.CRC.IBS$dataset <- "CRC"

## get top interactions which includes all the above shared taxa
## if number of assoc is > 10, pick top 10 assoc., else keep all assoc
gene_taxa_assoc <- CRC.gene.taxa.for.CRC.IBS
gene_taxa_assoc.filt <- data.frame()
if(dim(gene_taxa_assoc)[1] > 10){
  ## get top interactions which includes all the above shared taxa
   while(dim(gene_taxa_assoc.filt)[1] < 10 & dim(gene_taxa_assoc)[1] > 0){
    gene_taxa_assoc.dt <- setDT(gene_taxa_assoc)[order(BH_global), .SD[1], by=taxa]
    gene_taxa_assoc.filt <- rbind(gene_taxa_assoc.filt,as.data.frame(gene_taxa_assoc.dt))
    dim(gene_taxa_assoc.filt)
    gene_taxa_assoc <- gene_taxa_assoc[ !(gene_taxa_assoc$gene %in% gene_taxa_assoc.filt$gene & 
                                            gene_taxa_assoc$taxa %in% gene_taxa_assoc.filt$taxa),]
    dim(gene_taxa_assoc) 
    
   }
  dim(gene_taxa_assoc.filt) 
  ## if greater than 10 rows, subset to pick top 10
  CRC.gene.taxa.for.CRC.IBS <- gene_taxa_assoc.filt
  
}
dim(CRC.gene.taxa.for.CRC.IBS) 


IBS.gene.taxa.for.CRC.IBS <- IBS.lasso.stabsel[(IBS.lasso.stabsel$gene %in% IBS.genes.in.CRC.IBS.taxa.overlap & 
                                                  IBS.lasso.stabsel$taxa %in% CRC.IBS.taxa.overlap),]
dim(IBS.gene.taxa.for.CRC.IBS) 
IBS.gene.taxa.for.CRC.IBS$dataset <- "IBS"
## get top interactions which includes all the above shared taxa
## if number of assoc is > 10, pick top 10 assoc., else keep all assoc
gene_taxa_assoc <- IBS.gene.taxa.for.CRC.IBS
gene_taxa_assoc.filt <- data.frame()
if(dim(gene_taxa_assoc)[1] > 10){
  ## get top interactions which includes all the above shared taxa
   while(dim(gene_taxa_assoc.filt)[1] < 10 & dim(gene_taxa_assoc)[1] > 0){
    gene_taxa_assoc.dt <- setDT(gene_taxa_assoc)[order(BH_global), .SD[1], by=taxa]
    gene_taxa_assoc.filt <- rbind(gene_taxa_assoc.filt,as.data.frame(gene_taxa_assoc.dt))
    dim(gene_taxa_assoc.filt)
    gene_taxa_assoc <- gene_taxa_assoc[ !(gene_taxa_assoc$gene %in% gene_taxa_assoc.filt$gene & 
                                            gene_taxa_assoc$taxa %in% gene_taxa_assoc.filt$taxa),]
    dim(gene_taxa_assoc) 
    
   }
  dim(gene_taxa_assoc.filt) 
  ## if greater than 10 rows, subset to pick top 10
  IBS.gene.taxa.for.CRC.IBS <- gene_taxa_assoc.filt
  
}
dim(IBS.gene.taxa.for.CRC.IBS) 

## combine these assoc. for CRC ^ IBD 
combined_interactions <- rbind(CRC.gene.taxa.for.CRC.IBS,
                               IBS.gene.taxa.for.CRC.IBS
)
dim(combined_interactions) 

combined_interactions$taxa_name <- sapply(as.character(combined_interactions$taxa), pick_taxaname)
combined_interactions$shape_gene <- rep("gene",length(combined_interactions$gene))
combined_interactions$shape_taxa <- rep("taxa",length(combined_interactions$taxa))


combined_interactions$thickness <- abs(combined_interactions$spear_rho)
combined_interactions$sign <- sign(combined_interactions$ci.upper) 

## subset to relevant columns for viz
combined_interactions <- combined_interactions[c(1,2,13:20)]

################# Identify interactions for (IBD ^ IBS) taxa and plot them #############

## IBD ^ IBS
length(IBD.IBS.taxa.overlap) 
IBD.gene.taxa.for.IBD.IBS <-IBD.lasso.stabsel[(IBD.lasso.stabsel$gene %in% IBD.genes.in.IBD.IBS.taxa.overlap & 
                                                 IBD.lasso.stabsel$taxa %in% IBD.IBS.taxa.overlap),]
dim(IBD.gene.taxa.for.IBD.IBS) 
IBD.gene.taxa.for.IBD.IBS$dataset <- "IBD"

## get top interactions which includes all the above shared taxa
## if number of assoc is > 10, pick top 10 assoc., else keep all assoc
gene_taxa_assoc <- IBD.gene.taxa.for.IBD.IBS
gene_taxa_assoc.filt <- data.frame()
if(dim(gene_taxa_assoc)[1] > 10){
  ## get top interactions which includes all the above shared taxa
  while(dim(gene_taxa_assoc.filt)[1] < 10 & dim(gene_taxa_assoc)[1] > 0){
    gene_taxa_assoc.dt <- setDT(gene_taxa_assoc)[order(BH_global), .SD[1], by=taxa]
    gene_taxa_assoc.filt <- rbind(gene_taxa_assoc.filt,as.data.frame(gene_taxa_assoc.dt))
    dim(gene_taxa_assoc.filt)
    gene_taxa_assoc <- gene_taxa_assoc[ !(gene_taxa_assoc$gene %in% gene_taxa_assoc.filt$gene & 
                                            gene_taxa_assoc$taxa %in% gene_taxa_assoc.filt$taxa),]
    dim(gene_taxa_assoc) 
    
  }
  dim(gene_taxa_assoc.filt) 
  ## if greater than 10 rows, subset to pick top 10
  IBD.gene.taxa.for.IBD.IBS <- gene_taxa_assoc.filt
  
}
dim(IBD.gene.taxa.for.IBD.IBS) 

## Get IBS assoc. for IBD ^ IBS
IBS.gene.taxa.for.IBD.IBS <-IBS.lasso.stabsel[(IBS.lasso.stabsel$gene %in% IBS.genes.in.IBD.IBS.taxa.overlap & 
                                                 IBS.lasso.stabsel$taxa %in% IBD.IBS.taxa.overlap),]
IBS.gene.taxa.for.IBD.IBS$dataset <- "IBS"
dim(IBS.gene.taxa.for.IBD.IBS) 
## get top interactions which includes all the above shared taxa
## if number of assoc is > 10, pick top 10 assoc., else keep all assoc
gene_taxa_assoc <- IBS.gene.taxa.for.IBD.IBS
gene_taxa_assoc.filt <- data.frame()
if(dim(gene_taxa_assoc)[1] > 10){
  ## get top interactions which includes all the above shared taxa
   while(dim(gene_taxa_assoc.filt)[1] < 10 & dim(gene_taxa_assoc)[1] > 0){
    gene_taxa_assoc.dt <- setDT(gene_taxa_assoc)[order(BH_global), .SD[1], by=taxa]
    gene_taxa_assoc.filt <- rbind(gene_taxa_assoc.filt,as.data.frame(gene_taxa_assoc.dt))
    dim(gene_taxa_assoc.filt)
    gene_taxa_assoc <- gene_taxa_assoc[ !(gene_taxa_assoc$gene %in% gene_taxa_assoc.filt$gene & 
                                            gene_taxa_assoc$taxa %in% gene_taxa_assoc.filt$taxa),]
    dim(gene_taxa_assoc) 
    
   }
  dim(gene_taxa_assoc.filt)
  ## if greater than 10 rows, subset to pick top 10
  IBS.gene.taxa.for.IBD.IBS <- gene_taxa_assoc.filt
  
}
dim(IBS.gene.taxa.for.IBD.IBS) 


## combine these assoc. for CRC ^ IBD 
combined_interactions <- rbind(IBD.gene.taxa.for.IBD.IBS,
                               IBS.gene.taxa.for.IBD.IBS
)
dim(combined_interactions) 

combined_interactions$taxa_name <- sapply(as.character(combined_interactions$taxa), pick_taxaname)
combined_interactions$shape_gene <- rep("gene",length(combined_interactions$gene))
combined_interactions$shape_taxa <- rep("taxa",length(combined_interactions$taxa))

combined_interactions$thickness <- abs(combined_interactions$spear_rho)
combined_interactions$sign <- sign(combined_interactions$ci.upper) 

## subset to relevant columns for viz
combined_interactions <- combined_interactions[c(1,2,13:20)]

################# Identify interactions for (CRC ^ IBD ^ IBS) genes and plot them ##############  
length(CRC.IBD.IBS.genes)

## Get CRC assoc for these genes
CRC.gene.taxa.for.CRC.IBS.IBD.genes <- CRC.lasso.stabsel[(CRC.lasso.stabsel$gene %in% CRC.IBD.IBS.genes),]

CRC.gene.taxa.for.CRC.IBS.IBD.genes$dataset <- "CRC"
dim(CRC.gene.taxa.for.CRC.IBS.IBD.genes)


## Get IBD assoc for these genes
IBD.gene.taxa.for.CRC.IBS.IBD.genes <- IBD.lasso.stabsel[(IBD.lasso.stabsel$gene %in% CRC.IBD.IBS.genes),]

IBD.gene.taxa.for.CRC.IBS.IBD.genes$dataset <- "IBD"
dim(IBD.gene.taxa.for.CRC.IBS.IBD.genes)

## Get IBS assoc for these genes
IBS.gene.taxa.for.CRC.IBS.IBD.genes <- IBS.lasso.stabsel[(IBS.lasso.stabsel$gene %in% CRC.IBD.IBS.genes),]

IBS.gene.taxa.for.CRC.IBS.IBD.genes$dataset <- "IBS"
dim(IBS.gene.taxa.for.CRC.IBS.IBD.genes)

## combine these assoc. for CRC ^ IBD ^ IBS taxa for Viz.
combined_interactions <- rbind(CRC.gene.taxa.for.CRC.IBS.IBD.genes,
                               IBD.gene.taxa.for.CRC.IBS.IBD.genes,
                               IBS.gene.taxa.for.CRC.IBS.IBD.genes
)
dim(combined_interactions) 

combined_interactions$taxa_name <- sapply(as.character(combined_interactions$taxa), pick_taxaname)
combined_interactions$shape_gene <- rep("gene",length(combined_interactions$gene))
combined_interactions$shape_taxa <- rep("taxa",length(combined_interactions$taxa))

combined_interactions$thickness <- abs(combined_interactions$spear_rho)
combined_interactions$sign <- sign(combined_interactions$ci.upper) 

## subset to relevant columns for viz
combined_interactions <- combined_interactions[c(1,2,13:20)]

################# Identify interactions for (CRC ^ IBD) genes and plot them ##############  
length(CRC.IBD.genes.overlap) 


## Get CRC assoc for these genes
CRC.gene.taxa.for.CRC.IBD.genes <- CRC.lasso.stabsel[(CRC.lasso.stabsel$gene %in% CRC.IBD.genes.overlap),]

CRC.gene.taxa.for.CRC.IBD.genes$dataset <- "CRC"
dim(CRC.gene.taxa.for.CRC.IBD.genes)

## Get IBD assoc for these genes
IBD.gene.taxa.for.CRC.IBD.genes <- IBD.lasso.stabsel[(IBD.lasso.stabsel$gene %in% CRC.IBD.genes.overlap),]
IBD.gene.taxa.for.CRC.IBD.genes$dataset <- "IBD"
dim(IBD.gene.taxa.for.CRC.IBD.genes) 
IBD.gene.taxa.for.CRC.IBD.genes$gene[which(duplicated(IBD.gene.taxa.for.CRC.IBD.genes$gene))]


CRC.IBD.gene.taxa <- merge(CRC.gene.taxa.for.CRC.IBD.genes,IBD.gene.taxa.for.CRC.IBD.genes, by = c("gene"))
dim(CRC.IBD.gene.taxa)
View(CRC.IBD.gene.taxa)

## sort by CRC's BH_global, then IBD's BH_global
CRC.IBD.gene.taxa <- CRC.IBD.gene.taxa[order(CRC.IBD.gene.taxa$BH_global.x,CRC.IBD.gene.taxa$BH_global.y),]

## Let's pick top 10 unique genes from this table and 
## select assoc. from both datasets using these genes.
top10.genes <- unique(CRC.IBD.gene.taxa$gene)[1:10]

## Pick assoc. for these 10 genes
CRC.gene.taxa.for.CRC.IBD.genes <- CRC.gene.taxa.for.CRC.IBD.genes[(CRC.gene.taxa.for.CRC.IBD.genes$gene %in% top10.genes),]
dim(CRC.gene.taxa.for.CRC.IBD.genes)

## again, most genes are not duplicated in assoc. So let's take top 10
IBD.gene.taxa.for.CRC.IBD.genes <- IBD.gene.taxa.for.CRC.IBD.genes[(IBD.gene.taxa.for.CRC.IBD.genes$gene %in% top10.genes),]
dim(IBD.gene.taxa.for.CRC.IBD.genes)

## combine these assoc. for CRC ^ IBD genes for Viz.
combined_interactions <- rbind(CRC.gene.taxa.for.CRC.IBD.genes,
                               IBD.gene.taxa.for.CRC.IBD.genes
)
dim(combined_interactions)

combined_interactions$taxa_name <- sapply(as.character(combined_interactions$taxa), pick_taxaname)
combined_interactions$shape_gene <- rep("gene",length(combined_interactions$gene))
combined_interactions$shape_taxa <- rep("taxa",length(combined_interactions$taxa))


combined_interactions$thickness <- abs(combined_interactions$spear_rho)
combined_interactions$sign <- sign(combined_interactions$ci.upper)

## subset to relevant columns for viz
colnames(combined_interactions)
combined_interactions <- combined_interactions[c(1,2,13:20)]

################# Identify interactions for (CRC ^ IBS) genes and plot them ##############  
length(CRC.IBS.genes.overlap)

## Get CRC assoc for these genes
CRC.gene.taxa.for.CRC.IBS.genes <- CRC.lasso.stabsel[(CRC.lasso.stabsel$gene %in% CRC.IBS.genes.overlap),]

CRC.gene.taxa.for.CRC.IBS.genes$dataset <- "CRC"
dim(CRC.gene.taxa.for.CRC.IBS.genes)

## Get IBS assoc for these genes
IBS.gene.taxa.for.CRC.IBS.genes <- IBS.lasso.stabsel[(IBS.lasso.stabsel$gene %in% CRC.IBS.genes.overlap),]
IBS.gene.taxa.for.CRC.IBS.genes$dataset <- "IBS"
dim(IBS.gene.taxa.for.CRC.IBS.genes)
IBS.gene.taxa.for.CRC.IBS.genes$gene[which(duplicated(IBS.gene.taxa.for.CRC.IBS.genes$gene))]

CRC.IBS.gene.taxa <- merge(CRC.gene.taxa.for.CRC.IBS.genes,IBS.gene.taxa.for.CRC.IBS.genes, by = c("gene"))
dim(CRC.IBS.gene.taxa)
View(CRC.IBS.gene.taxa)

## sort by CRC's BH_global, then IBS's BH_global
CRC.IBS.gene.taxa <- CRC.IBS.gene.taxa[order(CRC.IBS.gene.taxa$BH_global.x,CRC.IBS.gene.taxa$BH_global.y),]

## Let's pick top 10 unique genes from this table and 
## select assoc. from both datasets using these genes.
top10.genes <- unique(CRC.IBS.gene.taxa$gene)[1:10]

## Pick assoc. for these 10 genes
CRC.gene.taxa.for.CRC.IBS.genes <- CRC.gene.taxa.for.CRC.IBS.genes[(CRC.gene.taxa.for.CRC.IBS.genes$gene %in% top10.genes),]
dim(CRC.gene.taxa.for.CRC.IBS.genes)

IBS.gene.taxa.for.CRC.IBS.genes <- IBS.gene.taxa.for.CRC.IBS.genes[(IBS.gene.taxa.for.CRC.IBS.genes$gene %in% top10.genes),]
dim(IBS.gene.taxa.for.CRC.IBS.genes)

## combine these assoc. for CRC ^ IBD genes for Viz.
combined_interactions <- rbind(CRC.gene.taxa.for.CRC.IBS.genes,
                               IBS.gene.taxa.for.CRC.IBS.genes
)
dim(combined_interactions)

combined_interactions$taxa_name <- sapply(as.character(combined_interactions$taxa), pick_taxaname)
combined_interactions$shape_gene <- rep("gene",length(combined_interactions$gene))
combined_interactions$shape_taxa <- rep("taxa",length(combined_interactions$taxa))

combined_interactions$thickness <- abs(combined_interactions$spear_rho)
combined_interactions$sign <- sign(combined_interactions$ci.upper)

## subset to relevant columns for viz
colnames(combined_interactions)
combined_interactions <- combined_interactions[c(1,2,13:20)]

################# Identify interactions for (IBD ^ IBS) genes and plot them ##############  
length(IBD.IBS.genes.overlap) 

## Get CRC assoc for these genes
IBD.gene.taxa.for.IBD.IBS.genes <- IBD.lasso.stabsel[(IBD.lasso.stabsel$gene %in% IBD.IBS.genes.overlap),]

IBD.gene.taxa.for.IBD.IBS.genes$dataset <- "IBD"
dim(IBD.gene.taxa.for.IBD.IBS.genes) 

IBD.IBS.gene.taxa <- merge(IBD.gene.taxa.for.IBD.IBS.genes,IBS.gene.taxa.for.IBD.IBS.genes, by = c("gene"))
dim(IBD.IBS.gene.taxa)
View(IBD.IBS.gene.taxa)

## sort by IBD's BH_global, then IBS's BH_global
IBD.IBS.gene.taxa <- IBD.IBS.gene.taxa[order(IBD.IBS.gene.taxa$BH_global.x,IBD.IBS.gene.taxa$BH_global.y),]

## Let's pick top 10 unique genes from this table and 
## select assoc. from both datasets using these genes.
top10.genes <- unique(IBD.IBS.gene.taxa$gene)[1:10]
top10.genes

## Pick assoc. for these 10 genes
IBD.gene.taxa.for.IBD.IBS.genes <- IBD.gene.taxa.for.IBD.IBS.genes[(IBD.gene.taxa.for.IBD.IBS.genes$gene %in% top10.genes),]
dim(IBD.gene.taxa.for.IBD.IBS.genes) 

IBS.gene.taxa.for.IBD.IBS.genes <- IBS.gene.taxa.for.IBD.IBS.genes[(IBS.gene.taxa.for.IBD.IBS.genes$gene %in% top10.genes),]
dim(IBS.gene.taxa.for.IBD.IBS.genes)

## combine these assoc. for CRC ^ IBD genes for Viz.
combined_interactions <- rbind(IBD.gene.taxa.for.IBD.IBS.genes,
                               IBS.gene.taxa.for.IBD.IBS.genes
)
dim(combined_interactions)

combined_interactions$taxa_name <- sapply(as.character(combined_interactions$taxa), pick_taxaname)
combined_interactions$shape_gene <- rep("gene",length(combined_interactions$gene))
combined_interactions$shape_taxa <- rep("taxa",length(combined_interactions$taxa))

combined_interactions$thickness <- abs(combined_interactions$spear_rho)
combined_interactions$sign <- sign(combined_interactions$ci.upper)

## subset to relevant columns for viz
colnames(combined_interactions)
combined_interactions <- combined_interactions[c(1,2,13:20)]

