## Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu

## Compute overlap between sparseCCA components across the diseases

## Initial setup 
rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
setwd("/path/to/working_directory")

## import libraries
library(msigdbr)
library(data.table)

################### Functions #####################
msigdb_enrichment <- function(pathway_DB, pathways, background_genes, genes_of_interest){
  
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
normalize_taxaname <- function(taxa){
  microbe_name <- taxa
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  return(microbe_name)
}

################### Set directories to read gene-taxa components for each dataset ############

msigdb_collection_name <- "C2_CP" 
dataset <- "CRC/case"
input_enrichment <- "sig_gene_taxa_components_padj"
CRC_comp_dir <- paste0("./data/analysis/",dataset,"/output_sparseCCA/grid_search/gene_taxa_components/",input_enrichment)

msigdb_collection_name <- "C2_CP" 
dataset <- "IBD/case"
input_enrichment <- "sig_gene_taxa_components_padj"
IBD_comp_dir <- paste0("./data/analysis/",dataset,"/output_sparseCCA/grid_search/gene_taxa_components/",input_enrichment)

msigdb_collection_name <- "C2_CP" 
dataset <- "IBS/case"
input_enrichment <- "sig_gene_taxa_components_padj"
IBS_comp_dir <- paste0("./data/analysis/",dataset,"/output_sparseCCA/grid_search/gene_taxa_components/",input_enrichment)

################## Compute gene set overlap across dieases ##################
CRC_comp_files <- list.files(CRC_comp_dir)
IBD_comp_files <- list.files(IBD_comp_dir)
IBS_comp_files <- list.files(IBS_comp_dir)
overlap_list <- list()
counter <- 0
for(i in 1:length(CRC_comp_files)){
  for(j in 1:length(IBD_comp_files)){
    for(k in 1:length(IBS_comp_files)){
      
      counter <- counter + 1
      CRC_genes <- read.table(paste0(CRC_comp_dir,"/", CRC_comp_files[i]),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
      CRC_genes <- CRC_genes$gene
      IBD_genes <- read.table(paste0(IBD_comp_dir,"/", IBD_comp_files[j]),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
      IBD_genes <- IBD_genes$gene
      IBS_genes <- read.table(paste0(IBS_comp_dir,"/", IBS_comp_files[k]),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
      IBS_genes <- IBS_genes$gene
      # overlap_genes <- Reduce(intersect, list(CRC_genes,IBD_genes,IBS_genes))
      ## Compute overlap using jaccrad index, i.e. intersection/union. Generalized for multiple sets.
      intersect_genes <- Reduce(intersect, list(CRC_genes,IBD_genes,IBS_genes))
      union_genes <- Reduce(union, list(CRC_genes,IBD_genes,IBS_genes))
      jaccard_index <- length(intersect_genes)/length(union_genes)
      
      overlap_df <- data.frame( CRC_comp = CRC_comp_files[i], IBD_comp = IBD_comp_files[j], IBS_comp = IBS_comp_files[k], 
                                count_genes =  length(intersect_genes), jaccard_index = jaccard_index, overlap_genes = paste(intersect_genes, collapse = ","))
      
      overlap_list[[counter]] <- overlap_df  
    }
  }
  
}

overlap_combined_df <- do.call(rbind.data.frame,overlap_list)
## sort by overlap_ratio
overlap_combined_df <- overlap_combined_df[order(overlap_combined_df$jaccard_index, decreasing = T),]

## list overlapping genes
cat(unlist(strsplit(as.character(overlap_combined_df$overlap_genes[1]), split=",")), sep = "\n")

################## Perform enrichment for overlapping gene sets using msigdb ##############

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

## Pathway DB of interest: KEGG and PID
path_DB <- c("CP:KEGG","CP:PID")
msigdb_C2_CP <- msigdb_C2_CP[(msigdb_C2_CP$gs_subcat %in% path_DB), ]
table(msigdb_C2_CP$gs_subcat)

## unique pathways in DBs of interest
msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
length(msigdb_C2_CP_unique_pathways) 
#382

## Perform enrichment
pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways
background_genes <- background_genes
## Check different overlp gene sets
genes_of_interest <- unlist(strsplit(as.character(overlap_combined_df$overlap_genes[1]), split = ","))

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