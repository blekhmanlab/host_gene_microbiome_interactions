### Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu


## This script performs postprocessing of lasso output and stability selection.


## Initial setup
rm(list=ls()) ##check workspace to make sure no precomputed data before clearing
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(RColorBrewer)
library(gridExtra)
library(corrplot)
library(data.table)


setwd("/path/to/working/dir") ## project directory

################# Functions ###################
## massage taxa names
massage_taxaname <- function(taxa){
  # microbe_name <- colnames(gene_x_microbe)[4] ## for debugging
  microbe_name <- taxa
  ## remove k__ or p__ etc. from the string
  microbe_name <- gsub("[A-z]__","", microbe_name)
  ## remove trailing ;
  microbe_name <- gsub("\\;+$","", microbe_name)
  return(microbe_name)
}


## Function to combine all the gene outputs together
get.combined.output <- function(filepath){
  filenames <- list.files(filepath, pattern="*.txt")
  summary_list <- list()
  count <- 1 
  for( file in filenames){
    # df <- read.table(file = paste0(filepath,"/",file), sep="\t",head=T, row.names = 1, check.names = F)
    df <- data.frame(fread(paste0(filepath,"/",file), sep="\t", head=T), row.names = 1, check.names = F)
    summary_list[[count]] <- df
    count <- count + 1
  }
  all_genes_result_df <- do.call(rbind,summary_list)
}

## correct for multiple hypothesis testing
MHT.correction <- function(df){
  ## Bonferonni
  assoc_Bonferonni <- p.adjust(df$pval, method = "bonferroni")
  df$Bonferonni_global <- assoc_Bonferonni
  ## BY
  assoc_BY <- p.adjust(df$pval, method="BY")
  df$BY_global <- assoc_BY
  ## BH 
  assoc_BH <- p.adjust(df$pval, method="BH")
  df$BH_global <- assoc_BH
  
  return(df)
}

create_scatterplot_IBS <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  ## add IBS_subytpe info before plotting
  IBS_subtype <- taxa.table$IBS_subtype + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,IBS_subtype,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    ## Note, colors are set by IBS_subtype
    geom_point(fill= col.fill[IBS_subtype], colour=border.col[IBS_subtype], shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  
  return(scatterplot)
}

create_scatterplot_Healthy <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}


## scatterplot function
create.scatterplot <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}

create.scatterplot.2 <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex)
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  condition <- taxa.table$condition + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,condition,pval,padj)
 
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill[condition], colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}

############### Combine MSI output files (SKIP IF ALREADY DONE) ################
## Input dataset
# dataset <- "CRC/case/"
# dataset <- "IBD/case/"
# dataset <- "IBS/case/"

# dataset <- "CRC/control/"
# dataset <- "IBD/control/"
dataset <- "IBS/control/"

combined.output <- get.combined.output(paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V5/output_MSI/"))
dim(combined.output)

# write.table(combined.output, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V5/combined_output/combined_lasso_output.txt"), sep="\t", row.names = T, col.names = NA )

############## Input combined output ###############
# dataset <- "CRC/case/"
# dataset <- "IBD/case/"
# dataset <- "IBS/case/"

# dataset <- "CRC/control/"
# dataset <- "IBD/control/"
dataset <- "IBS/control/"

combined.output.filepath <- paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V5/combined_output/combined_lasso_output.txt")
## read.table() takes too long! Replace with fread(), see below.
## fread() to read large file. 
combined.output <- data.frame(fread(combined.output.filepath, sep="\t", head=T), row.names = 1)
dim(combined.output)

########### Perform MHT correction on p-values for all assoc. ############

combined.output <- MHT.correction(combined.output)
## Order output by BH_global
combined.output <- combined.output[order(combined.output$BH_global),]

## write to file
# filename <- "combined_lasso_output_MHT_corrected.txt"
# write.table(combined.output, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V5/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

################# SKIP IF DONE -- Combine output from MSI stability selection run ##############
## Combine all the gene outputs together
stabsel_path <- "stability_selection_V5"
stability.combined.output <- get.combined.output(paste0("./data/analysis/",dataset,"/",stabsel_path,"/output_MSI/"))
dim(stability.combined.output)

length(unique(stability.combined.output$gene_name)) 

## Number of associations with a stability selected covariate (taxa or gender or condition)
dim(stability.combined.output[stability.combined.output$taxa.selected != "None",])

## Number of genes with stability selected covariate (taxa or gender or condition)
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None",]$gene_name))

## Number of genes with stability selected taxa
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None" & stability.combined.output$taxa.selected != "gender",]$gene_name)) 

## write to file
# stabsel.combined.filename <- "combined_stabsel.txt"
# write.table(stability.combined.output, file = paste0("./data/analysis/",dataset,"/",stabsel_path,"/combined_output/", stabsel.combined.filename), sep="\t", row.names = T, col.names = NA )

################# Test overlap with stability selected output ###############

## Input dataset and stabsel output filename
stabsel_path <- "stability_selection_V5"
stabsel_file <- "combined_stabsel.txt"

stabsel_output_filepath <- paste0("./data/analysis/",dataset,"/",stabsel_path,"/combined_output/",stabsel_file)
stabsel.output <- read.table(stabsel_output_filepath,sep="\t",head=T, row.names = 1, check.names = F)
dim(stabsel.output)

## massage column names of stability selected output for downstream processing
colnames(stabsel.output) <- c("gene","taxa")

## Identify stability selected associations
stabsel.output.filt <- stabsel.output[stabsel.output$taxa != "None",] 
dim(stabsel.output.filt)

## Identify unique stability selected genes
genes.stabsel <- unique(stabsel.output.filt$gene)
length(genes.stabsel)

## Identify unique stability selected predictors
taxa.stabsel <- unique(stabsel.output.filt$taxa)
length(taxa.stabsel)

## Overlap between combined.output.filt and stability selected assoc
overlap.assoc.stabsel <- merge(combined.output,stabsel.output.filt, by = c("gene","taxa"))
overlap.assoc.stabsel <- overlap.assoc.stabsel[order(overlap.assoc.stabsel$BH_global),]
dim(overlap.assoc.stabsel)

length(unique(overlap.assoc.stabsel$gene))

length(unique(overlap.assoc.stabsel$taxa))

## write to file before filtering by FDR cutoffs
# filename <- "combined_output_MHT_stabsel.txt"
# write.table(overlap.assoc.stabsel, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V5/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

## How many stability selected assoc. also significant at FDR < 0.1
overlap.assoc.stabsel.BH.0.1 <- overlap.assoc.stabsel[overlap.assoc.stabsel$BH_global < 0.1,]
dim(overlap.assoc.stabsel.BH.0.1)

length(unique(overlap.assoc.stabsel.BH.0.1$gene))

length(unique(overlap.assoc.stabsel.BH.0.1$taxa)) 

############# Filter out non-taxa associations (e.g. gene-gender assoc.)######################

combined.output <- overlap.assoc.stabsel.BH.0.1

## Identify gender associations
select <- which(combined.output$taxa =="gender")
## Spit out gender associated output to a file for later inspection
gender.assoc <- combined.output[select,]; dim(gender.assoc)

length(unique(gender.assoc$gene))

## write to file
# write.table(gender.assoc, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V5/combined_output/gene_gender_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )

## filter out gender associations.
combined.output.filt <- combined.output[-select,]
dim(combined.output.filt)

## Additional gene-covariate filtering for IBD or IBS
if(dataset == "IBD/case/" || dataset == "IBS/case/"){
  ## set subtype label
  subtype_covariate <- "IBS_subtype" #IBD_subtype
  select <- which(combined.output.filt$taxa == subtype_covariate)
  condition.assoc <- combined.output.filt[select,]; dim(condition.assoc)
  
  length(unique(condition.assoc$gene)) 
  
  # write gene-condition associations to file
  # write.table(condition.assoc, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V2/combined_output/gene_IBS_subtype_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )
  
  # #keep non-condition associations.
  if(length(select) != 0)
    combined.output.filt <- combined.output.filt[-select,]; dim(combined.output.filt)
}


## no. of unique genes
length(unique(combined.output.filt$gene))

## no. of unique taxa
length(unique(combined.output.filt$taxa))

## set datatset and filename and write to file
# dataset <- "IBS/control/"
# filename <- "IBS_control_gene_taxa_FDR_0.1.txt"
# write.table(combined.output.filt, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V5/combined_output/",filename), sep="\t", row.names = T, col.names = NA )


