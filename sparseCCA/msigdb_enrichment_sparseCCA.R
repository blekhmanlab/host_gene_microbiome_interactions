## Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu

## This script computes enrichment of pathways for genes in significant components of sparseCCA analysis
## This uses collection of pathways from MSigDB. 

## Perform this for each dataset's sparse CCA output (CRC, IBD,IBS)
## For each dataset, we perform erichment per component, correct for MHT (BH-FDR), pool pathways across components,
## sort by FDR, and filter duplicated pathways.This would result in assignment of smallest FDR per pathway.
## We repeat this for case and control per dataset. 
## We use CompGO's odd-ratio based differential enrichment strategy to compute case-specific pathways. 


## Initial setup 
rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
setwd("/path/to/working_directory")

# install.packages("msigdbr")
library(msigdbr)
library(data.table)

############### Functions ################

## Inputs:
## A database of pathways
## Background genes
## Genes list of interest

## Output: 
## Enrichment result

msigdb_enrichment <- function(pathway_DB, pathways, background_genes, genes_of_interest){
  
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
    
    ## m, #overlapping genes in this pathway and background genes for CRC
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
    
    ## Extract list of genes in the genes of interest that are included in the pathway. 
    gene_names_in_pathway = noquote(paste(intersect(pathway_gene_set,genes_of_interest), collapse = ","))
   
     ## k, total #genes in genes list of interest
    k <- length(genes_of_interest)
    k
    
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

## Do enrichment per component
perform_enrichment <- function(input_dir, filenames, msigdb_collection_name, pathway_DB, pathways, background_genes, output_dir){
  enriched_list <- list()
  count <- 0
  
  ## perform enrichment for each component 
  for(i in filenames){
    
    ## debug
    # i <- "gene_taxa_component_1.txt"
    
    count <- count + 1
    
    print(paste0("Enrichment for component: ", i));flush.console()
    ## load genes list for a given component
    sparseCCA_genes <- read.table(paste0(input_dir,"/",i),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
    sparseCCA_genes <- sparseCCA_genes$gene
    genes_of_interest <- sparseCCA_genes
    
    ## perform enrichment using pathways in msigdb collection
    msigdb_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
    ## convert the list to dataframe
    msigdb_pathways_df <- do.call(rbind.data.frame, msigdb_pathways)
    ## Add current component as a column to enrichment result
    msigdb_pathways_df$Component <- rep(paste0(i),nrow(msigdb_pathways_df))
    ## Add collection name as a column to the enrichemnt result
    msigdb_pathways_df$Collection <- msigdb_collection_name
    
    ## Make last column(i.e., component) as first column
    msigdb_pathways_df <- msigdb_pathways_df[,c(ncol(msigdb_pathways_df)-1, ncol(msigdb_pathways_df),1:(ncol(msigdb_pathways_df)-2))]
    
    ## update column name to reflect component, except pathway colname
    # colnames(msigdb_pathways_df[,-3]) <- paste(colnames(msigdb_pathways_df[,-3]), count, sep = "_")
    
    ## sort pathways by pval
    msigdb_pathways_df <- msigdb_pathways_df[order(msigdb_pathways_df$p_val, decreasing = F),]
    length(msigdb_pathways_df$pathway)
    
    ## save dataframe of pathways for a component to file
    if(!is.null(output_dir)){
      filename <- strsplit(i,"\\.")[[1]][1]
      write.table(msigdb_pathways_df, file = paste0(output_dir,"/msigdb_",filename,".txt") , sep="\t", row.names = F)
    }
    ## add df to list
    enriched_list[[count]] <- msigdb_pathways_df
  }
  
  ## This is the list of dataframes, where each dataframe holds enrichment result for a component.  
  return(enriched_list)
}

############### Input background genes (case and control) ###############
# Provide universe/background genes used as input for sparseCCA analysis

## case genes
filename <- "./data/clean/CRC/gene_expr/tumor_protein_coding_vsd_SDfilt_t.txt"
genes <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
genes <- t(genes)
dim(genes)
CRC.case.background.genes  <- rownames(genes)
CRC.case.background.genes <- CRC.case.background.genes[order(CRC.case.background.genes)]
length(CRC.case.background.genes)

# control genes
filename <- "./data/clean/CRC/gene_expr/normal_protein_coding_vsd_SDfilt_t.txt"
genes_ctrl <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
genes_ctrl <- t(genes_ctrl)
dim(genes_ctrl)
CRC.ctrl.background.genes  <- rownames(genes_ctrl)
CRC.ctrl.background.genes <- CRC.ctrl.background.genes[order(CRC.ctrl.background.genes)]
length(CRC.ctrl.background.genes)

############### Get MSigDB pathways #################

## Let's focus first on C2 subcollection canonical pathways 
## http://software.broadinstitute.org/gsea/msigdb/collection_details.jsp#C2
## This includes gene sets from BioCarta, KEGG, PID, and Reactome. 

msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
class(msigdb_C2) #[1] "tbl_df"     "tbl"        "data.frame"
msigdb_C2 <- as.data.frame(msigdb_C2)
table(msigdb_C2$gs_subcat)
# CGP          CP CP:BIOCARTA     CP:KEGG      CP:PID CP:REACTOME 
# 363007        3697        4775       12783        8058       85647
## Only keep canonical pathways, or in other words, remove CGP (Chemical and genetic perturbations)

msigdb_C2_CP <- msigdb_C2[msigdb_C2$gs_subcat!= "CGP",]
table(msigdb_C2_CP$gs_subcat)
# CP CP:BIOCARTA     CP:KEGG      CP:PID CP:REACTOME 
# 3697        4775       12783        8058       85647 
dim(msigdb_C2_CP)
# [1] 114960      9
length(unique(msigdb_C2_CP$gs_name))
# 2199 pathways

## Only keep pathway DBs that we want to test gene sets against 
table(msigdb_C2_CP$gs_subcat)
# CP CP:BIOCARTA     CP:KEGG      CP:PID CP:REACTOME 
# 3697        4775       12783        8058       85647 

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

############### Case: Set filtering criteria for all pathway DBs ##########

##(25,300,5)
dataset <- "CRC/case"
min_genes <- 25
max_genes <- 300
overlap_genes <- 5

############### Case enrichment using msigdb_C2_CP ############

## Perform enrichment per component. 
#Set msigdb_collection_name, input_dir, filenames, pathway_DB, pathways, background_genes, output_dir
msigdb_collection_name <- "C2_CP" 
dataset <- "CRC/case"
input_enrichment <- "sig_gene_taxa_components_0.2_0.15_padj"
input_dir <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/gene_taxa_components/",input_enrichment)
filenames <- list.files(input_dir)
stopifnot(length(filenames) > 0) # make sure the input dir path is set correctly
pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways
background_genes <- CRC.case.background.genes
output_dir <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/enrichment_sparseCCA_msigDB/",msigdb_collection_name,"/per_component_kegg_pid_reactome/before_MHT/25_300_5/")
## create the output dir if it doesnot exist
ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = T), FALSE)

## Perform enrichment
enrichment_list <- perform_enrichment(input_dir, filenames, 
                                      msigdb_collection_name, pathway_DB, pathways, background_genes, 
                                      output_dir)

## Pool enrichment output across components.
enrichment_case <- do.call(rbind,enrichment_list)
dim(enrichment_case)
## Any duplicates?
any(duplicated(enrichment_case$pathway))

## order by p-value
enrichment_case <- enrichment_case[order(enrichment_case$p_val),]

##filter duplicates (keeps pathways with smalles p-value across components)
enrichment_case_no_dups <- enrichment_case[!duplicated(enrichment_case$pathway),]
dim(enrichment_case_no_dups) 

## MHT correction -- FDR
enrichment_case_no_dups$p_adj <- p.adjust(enrichment_case_no_dups$p_val, method = "BH")
enrichment_case_no_dups <- enrichment_case_no_dups[order(enrichment_case_no_dups$p_adj),]
length(which(enrichment_case_no_dups$p_adj < 0.1))

case_kegg_pid <- enrichment_case_no_dups

## write to file
# filename <- paste0("CRC_case_msigdb.txt")
# filepath <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/final_gene_taxa_components/enrichment_sparseCCA_msigDB/",msigdb_collection_name,"/",filename)
# write.table(enrichment_case, file = paste0(filepath) , sep="\t", row.names = F)

############### CompGO case: log odds ratio based z-score ############
## Based on the idea of identifying differential pathways between experiments stated in the CompGO paper
## Waardenberg et al., BMC Bioinformatics, 2015
## https://link.springer.com/article/10.1186/s12859-015-0701-2

## Compute SE(log OR) = sqrt(1/n1 + 1/n2 + 1/n3 + 1/n4) per pathway
case_kegg_pid$SE_logOR <- sqrt(1/case_kegg_pid$cont_n1 + 1/case_kegg_pid$cont_n2 + 1/case_kegg_pid$cont_n3 + 1/case_kegg_pid$cont_n4)

## Compute z-score = log(OR)/SE(logOR)
case_kegg_pid$z_logOR <- log(case_kegg_pid$odds_ratio)/case_kegg_pid$SE_logOR

############### Control: Set filtering criteria for all pathway DBs ##########
##(25,300,5)
dataset <- "CRC/control"
min_genes <- 25
max_genes <- 300
overlap_genes <- 5


############### Control enrichment using msigdb_C2_CP ###########

#Set msigdb_collection_name, input_dir, filenames, pathway_DB, pathways, background_genes, output_dir
msigdb_collection_name <- "C2_CP" 
dataset <- "CRC/control"
input_dir <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/gene_taxa_components/sig_gene_taxa_components_0.266666666666667_0.344444444444444_padj")
filenames <- list.files(input_dir)
stopifnot(length(filenames) > 0) # make sure the input dir path is set correctly
pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways
background_genes <- CRC.ctrl.background.genes
output_dir <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/enrichment_sparseCCA_msigDB/",msigdb_collection_name,"/per_component_kegg_pid/before_MHT/25_300_5/")
## create the output dir if it doesnot exist
ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = T), FALSE)


## Perform enrichment
# output_dir <- NULL ## to prevent overwriting per component enrichment result
enrichment_list <- perform_enrichment(input_dir, filenames, 
                                      msigdb_collection_name, pathway_DB, pathways, background_genes, 
                                      output_dir)
enrichment_control <- do.call(rbind,enrichment_list)
dim(enrichment_control) 

## Any duplicates?
any(duplicated(enrichment_control$pathway)) #T

## order by p-value
enrichment_control <- enrichment_control[order(enrichment_control$p_val),]

##filter duplicates (keeps pathways with smalles p-value across components)
enrichment_control_no_dups <- enrichment_control[!duplicated(enrichment_control$pathway),]
dim(enrichment_control_no_dups) 

## MHT correction -- FDR
enrichment_control_no_dups$p_adj <- p.adjust(enrichment_control_no_dups$p_val, method = "BH")
enrichment_control_no_dups <- enrichment_control_no_dups[order(enrichment_control_no_dups$p_adj),]
length(which(enrichment_control_no_dups$p_adj < 0.1))

control_kegg_pid <- enrichment_control_no_dups


## write to file
# filename <- paste0("CRC_control_msigdb.txt")
# filepath <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/final_gene_taxa_components/enrichment_sparseCCA_msigDB/",msigdb_collection_name,"/",filename)
# write.table(enrichment_control, file = paste0(filepath) , sep="\t", row.names = F)

############### CompGO control: log odds ratio based z-score ############
## Based on the idea of identifying differential pathways between experiments stated in the CompGO paper
## Waardenberg et al., BMC Bioinformatics, 2015
## https://link.springer.com/article/10.1186/s12859-015-0701-2

## Note, the log odds ratio has an approximately normal distribution, hence we can compute the 
## standard error and confidence interval for log odds ratio. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1127651/)

## Compute SE(log OR) = sqrt(1/n1 + 1/n2 + 1/n3 + 1/n4) per pathway
## Apart from CompGO, formula for SE(logOR) is confirmed in Morris et al. 1988, 
## and Bland et al. 2000 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1127651/)
control_kegg_pid$SE_logOR <- sqrt(1/control_kegg_pid$cont_n1 + 1/control_kegg_pid$cont_n2 + 1/control_kegg_pid$cont_n3 + 1/control_kegg_pid$cont_n4)

## Compute z-score = log(OR)/SE(logOR)
control_kegg_pid$z_logOR <- log(control_kegg_pid$odds_ratio)/control_kegg_pid$SE_logOR

############### CompGO Case - Control: differential enrichment using comparative log odds ratio ##############
## Filter case and control pathways by FDR cut-off (0.1) and 
## perform differential enrichment using z-score method.

## Check how many case and control pathways overlap at FDR < 0.1
case_kegg_pid_FDR_0.1 <- case_kegg_pid[(case_kegg_pid$p_adj < 0.1),]
control_kegg_pid_FDR_0.1 <- control_kegg_pid[(control_kegg_pid$p_adj < 0.1),]

## Overlapping pathways between case and control filtered at FDR < 0.1
case_control_overlap_FDR_0.1 <- merge(case_kegg_pid_FDR_0.1,
                                      control_kegg_pid_FDR_0.1,
                                      by = "pathway")
dim(case_control_overlap_FDR_0.1) 

## How many case-only pathways at FDR < 0.1 (i.e. only found in case, not in control at FDR < 0.1)
case_only_path_FDR_0.1 <- case_kegg_pid_FDR_0.1[!(case_kegg_pid_FDR_0.1$pathway %in% control_kegg_pid_FDR_0.1$pathway),]
dim(case_only_path_FDR_0.1)[1]

## Compare overlapping pathways
## compute differential z-score = log(OR)_case - log(OR_ctrl)/ sqrt(SE_logOR_case^2 + SE_logOR_ctrl^2)
case_control_overlap_FDR_0.1$z_case_ctrl <- (log(case_control_overlap_FDR_0.1$odds_ratio.x) - log(case_control_overlap_FDR_0.1$odds_ratio.y))/sqrt((case_control_overlap_FDR_0.1$SE_logOR.x)^2 + (case_control_overlap_FDR_0.1$SE_logOR.y)^2)

## Assuming normality, compute p-value corresponding to the z-score
## Let's compute p-val for one-tailed test, i.e. to test if pathway is over-enriched in case vs control
## For 1-sided test, first only keep pathways that have +ve z-score
## since -ve z-score implies that the pathway could be potentially more enriched in control.
case_control_overlap_FDR_0.1_1_tail <- case_control_overlap_FDR_0.1[(case_control_overlap_FDR_0.1$z_case_ctrl > 0),]
dim(case_control_overlap_FDR_0.1_1_tail) #[1] 70 38
case_control_overlap_FDR_0.1_1_tail$z_pval <- pnorm(-abs(case_control_overlap_FDR_0.1_1_tail$z_case_ctrl))
length(which(case_control_overlap_FDR_0.1_1_tail$z_pval < 0.05)) #19
length(which(case_control_overlap_FDR_0.1_1_tail$z_pval < 0.1)) #26

##sort by pval
case_control_overlap_FDR_0.1_1_tail <- case_control_overlap_FDR_0.1_1_tail[order(case_control_overlap_FDR_0.1_1_tail$z_pval, decreasing = F),]

## FDR adjustment on pvalues
case_control_overlap_FDR_0.1_1_tail$z_FDR <- p.adjust(case_control_overlap_FDR_0.1_1_tail$z_pval, method = "BH")
length(which(case_control_overlap_FDR_0.1_1_tail$z_FDR < 0.2)) #19
length(which(case_control_overlap_FDR_0.1_1_tail$z_FDR < 0.2 & case_control_overlap_FDR_0.1_1_tail$z_pval < 0.05)) #19

dataset <- "CRC/case"
filename <- "CRC_case_control_overlap_CompGO_diff_pathway_FDR_0.1.txt" 
filepath <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/enrichment_sparseCCA_msigDB/C2_CP/case_minus_control_kegg_pid/FDR_cutoff/")
## create the output dir if it doesnot exist
ifelse(!dir.exists(filepath), dir.create(filepath, recursive = T), FALSE)
# write.table(case_control_overlap_FDR_0.1_1_tail, file = paste0(filepath,filename) , sep="\t", row.names = F)

##list DE pathways
DE_case_pathways <- as.character(case_control_overlap_FDR_0.1_1_tail[which(case_control_overlap_FDR_0.1_1_tail$z_FDR < 0.2 & case_control_overlap_FDR_0.1_1_tail$p_adj.x < 0.1),]$pathway)

## compute gene set overlap coeff and jaccard index for common pathways
case_control_overlap_FDR_0.1_1_tail$case_control_OC <- NA
case_control_overlap_FDR_0.1_1_tail$case_control_JC <- NA

for(i in 1:nrow(case_control_overlap_FDR_0.1_1_tail)){
  
  gene_set_case <-  case_control_overlap_FDR_0.1_1_tail[i,]$gene_names_in_pathway.x
  gene_set_case <- unlist(strsplit(gene_set_case, split = ","))
  
  gene_set_control <- as.character(case_control_overlap_FDR_0.1_1_tail[i,]$gene_names_in_pathway.y)
  gene_set_control <- unlist(strsplit(gene_set_control, split = ","))
  
  ## compute overlap coefficient between these two gene sets
  ## Overlap coefficient is defined as: the size of the intersection divided by the smaller of the size of the two sets.
  ## overlap(X,Y) = |X ^ Y|/min(|X|,|Y|)
  ## src: https://en.wikipedia.org/wiki/Overlap_coefficient
  case_control_intersect <- length(intersect(gene_set_case,gene_set_control)) 
  case_control_union <- length(union(gene_set_case,gene_set_control))
  overlap_coeff <- case_control_intersect/min(length(gene_set_case),length(gene_set_control))
  jaccard_coeff <- case_control_intersect/case_control_union
  case_control_overlap_FDR_0.1_1_tail$case_control_OC[i] <- overlap_coeff
  case_control_overlap_FDR_0.1_1_tail$case_control_JC[i] <- jaccard_coeff
  
}

## Select case pathwys that satisfy the following:
## 1. case-only pathways at FDR < 0.1, i.e. case pathways not found in control
## 2. pathways differentially enriched in case vs control at FDR < 0.2 and pathway enrichment score in case < 0.1 
select <- c(as.character(case_only_path_FDR_0.1$pathway), DE_case_pathways)
case_specific_pathways <- case_kegg_pid_FDR_0.1[(case_kegg_pid_FDR_0.1$pathway %in% select),]
dim(case_specific_pathways)

dataset <- "CRC/case"
filename <- "CRC_case_specific_CompGO_path_FDR_0.1_filter_25_300_5.txt"
filepath <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/enrichment_sparseCCA_msigDB/C2_CP/case_minus_control_kegg_pid/FDR_cutoff/")
## create the output dir if it doesnot exist
ifelse(!dir.exists(filepath), dir.create(filepath, recursive = T), FALSE)
# write.table(case_specific_pathways, file = paste0(filepath,filename) , sep="\t", row.names = F)
