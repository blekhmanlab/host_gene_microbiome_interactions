## Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu

## This script performs enrichment for disease-specific genes found significantly 
## associated with a taxa. 

## Initial setup 
rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
setwd("/path/to/working_directory")

############### Function for enrichment ##############
perform_enrichment<- function(genes,db){
  
  enriched_pathways <- enrichr(case.only.genes,db)
  ## convert the list of pathway vectors into dataframe
  enriched_pathways_df <- do.call(rbind.data.frame, enriched_pathways)
  
  ## Make current rowname (without the trailing .1,.2, etc) as a column called Database
  db_name <- rownames(enriched_pathways_df)
  db_name <- gsub("\\.\\d+","",db_name)
  enriched_pathways_df$Database <- db_name
  ## drop rownames 
  rownames(enriched_pathways_df) <- c()
  ## Make last database as first column
  enriched_pathways_df <- enriched_pathways_df[,c(ncol(enriched_pathways_df),1:(ncol(enriched_pathways_df)-1))]
  
  ## sort pathways by adjusted pval
  enriched_pathways_df <- enriched_pathways_df[order(enriched_pathways_df$Adjusted.P.value, decreasing = F),]
  
  enriched_pathways_df
}

msigdb_enrichment <- function(pathway_DB, pathways, background_genes, genes_of_interest){
  
  ## set min and max genes to allow pathways for testing
  min_genes <- 25
  max_genes <- 85
  overlap_genes <- 5
  
  enrichment_list <- list()
  for(i in 1:length(pathways)){
    
    ## debug
    # i <- 56
    # pathway_DB <- msigdb_C2_CP
    # pathways <- msigdb_C2_CP_unique_pathways
    # background_genes <- background_genes
    # ## for now, check 1st overlap set
    # genes_of_interest <- unlist(strsplit(as.character(overlap_combined_df$overlap_genes[1]), split = ","))
    # 
    ## debug
    
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
    ## Good explanation here: https://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
    ## Matches orientation of contingency rows and col of DAVID: https://david.ncifcrf.gov/content.jsp?file=functional_annotation.html
    contingency_table <- matrix(c(x, k-x, m-x, n-(k-x)), #matrix is filled along the column by default
                                nrow = 2, ncol = 2, 
                                dimnames = list(c("In_pathway","Not_in_pathway"),
                                                c("Genes_of_interest","Genes_NOT_of_interest"))
    )
    
    contingency_table
    #                       Genes_of_interest Genes_NOT_of_interest
    # In_pathway                     1                    10
    # Not_in_pathway               838                 11664
    
    fisher_result <- fisher.test( contingency_table, alternative = "greater")
    #print(fisher_result)
    # Fisher's Exact Test for Count Data
    # data:  contingency_table
    # p-value = 8.934e-06
    # alternative hypothesis: true odds ratio is greater than 1
    # 95 percent confidence interval:
    #   2.620139      Inf
    # sample estimates:
    #   odds ratio 
    # 4.555846 
    
    ## For debugging, check with hypergeometric test
    # phyper_pval <- phyper(x-1, m-x, n-(k-x), k, lower.tail = F)
    # phyper_pval
    # 1-sum(dhyper(0:(x-1),m-x, n-(k-x), k))
    
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

############### Input background genes #################
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

############### Input genes of interest (case gene-taxa interactions) #############
dataset <- "CRC/case"
filename <- "CRC_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("data/analysis/",dataset,"/output_lasso_hdi_loocv_V5/combined_output/",filename)
CRC.lasso.stabsel <- read.table(filepath,sep="\t",head=T, row.names = 1, check.names = F)
dim(CRC.lasso.stabsel)
CRC.case.genes <- unique(CRC.lasso.stabsel$gene)
length(CRC.case.genes)

dataset <- "IBD/case"
filename <- "IBD_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("data/analysis/",dataset,"/output_lasso_hdi_loocv_V5/combined_output/",filename)
IBD.lasso.stabsel <- read.table(filepath,sep="\t",head=T, row.names = 1, check.names = F)
dim(IBD.lasso.stabsel)
IBD.case.genes <- unique(IBD.lasso.stabsel$gene)
length(IBD.case.genes)

dataset <- "IBS/case"
filename <- "IBS_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("data/analysis/",dataset,"/output_lasso_hdi_loocv_V5/combined_output/",filename)
IBS.lasso.stabsel <- read.table(filepath,sep="\t",head=T, row.names = 1, check.names = F)
dim(IBS.lasso.stabsel) 
IBS.case.genes <- unique(IBS.lasso.stabsel$gene)
length(IBS.case.genes)

############### Input genes of interest (case-specific genes = case - control) #############
## Here we read case-specific genes

dataset <- "CRC/case"
filename <- "CRC_case_specific_genes_assoc_w_taxa.txt"
filepath <- paste0("data/analysis/",dataset,"/output_lasso_hdi_loocv_V5/combined_output/",filename)
CRC.case.genes <- scan(filepath, what="", sep="\n")
length(CRC.case.genes)

dataset <- "IBD/case"
filename <- "IBD_case_specific_genes_assoc_w_taxa.txt"
filepath <- paste0("data/analysis/",dataset,"/output_lasso_hdi_loocv_V5/combined_output/",filename)
IBD.case.genes <- scan(filepath, what="", sep="\n")
length(IBD.case.genes)

dataset <- "IBS/case"
filename <- "IBS_case_specific_genes_assoc_w_taxa.txt"
filepath <- paste0("data/analysis/",dataset,"/output_lasso_hdi_loocv_V5/combined_output/",filename)
IBS.case.genes <- scan(filepath, what="", sep="\n")
length(IBS.case.genes)

############### Get MSigDB pathways #################

## Let's focus first on C2 subcollection canonical pathways 
## http://software.broadinstitute.org/gsea/msigdb/collection_details.jsp#C2
## This includes gene sets from BioCarta, KEGG, PID, and Reactome. 

msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
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

## Pathway DB of interest
path_DB <- c("CP:KEGG","CP:PID", "CP:REACTOME")
msigdb_C2_CP <- msigdb_C2_CP[(msigdb_C2_CP$gs_subcat %in% path_DB), ]
table(msigdb_C2_CP$gs_subcat)
# CP:KEGG      CP:PID CP:REACTOME 
# 12783        8058       85647 

## unique pathways in DBs of interest
msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
length(msigdb_C2_CP_unique_pathways) 
#1881
############## CRC case enrichment ##################
## Set enrichment call parameters
pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways
background_genes <- CRC.case.background.genes
genes_of_interest <- CRC.case.genes

length(background_genes)
length(genes_of_interest)

##Call enrichment function 
CRC_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
## convert the list to dataframe
CRC_pathways <- do.call(rbind.data.frame, CRC_pathways)
dim(CRC_pathways) 

##sort by pval
CRC_pathways <- CRC_pathways[order(CRC_pathways$p_val),]
## MHT correction (FDR)
CRC_pathways$p_adj <- p.adjust(CRC_pathways$p_val, method = "BH")

length(which(CRC_pathways$p_adj < 0.2))

CRC_pathways <- CRC_pathways[CRC_pathways$p_adj < 0.2,]
dim(CRC_pathways) 

############## IBD case enrichment ##################
## Set enrichment call parameters
pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways
background_genes <- IBD.case.background.genes
genes_of_interest <- IBD.case.genes

length(background_genes)
length(genes_of_interest)

##Call enrichment function 
IBD_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
## convert the list to dataframe
IBD_pathways <- do.call(rbind.data.frame, IBD_pathways)
dim(IBD_pathways) 

##sort by pval
IBD_pathways <- IBD_pathways[order(IBD_pathways$p_val),]
## MHT correction (FDR)
IBD_pathways$p_adj <- p.adjust(IBD_pathways$p_val, method = "BH")

length(which(IBD_pathways$p_adj < 0.2)) 

IBD_pathways <- IBD_pathways[IBD_pathways$p_adj < 0.2,]
dim(IBD_pathways)

############## IBS case enrichment ##################
## Set enrichment call parameters
pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways
background_genes <- IBS.case.background.genes
genes_of_interest <- IBS.case.genes

length(background_genes)
length(genes_of_interest)

##Call enrichment function 
IBS_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
## convert the list to dataframe
IBS_pathways <- do.call(rbind.data.frame, IBS_pathways)
dim(IBS_pathways) 

##sort by pval
IBS_pathways <- IBS_pathways[order(IBS_pathways$p_val),]
## MHT correction (FDR)
IBS_pathways$p_adj <- p.adjust(IBS_pathways$p_val, method = "BH")

length(which(IBS_pathways$p_adj < 0.2)) 

IBS_pathways <- IBS_pathways[IBS_pathways$p_adj < 0.2,]
dim(IBS_pathways)