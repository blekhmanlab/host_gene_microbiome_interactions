### Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu

## This script performs stability selection for lasso model to identify
## robustly associated microbes with host genes.
## Here we use the "stabs" R package's implementation of 
## stability selection.
## Stabs paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0575-3
## Stabs manual: https://cran.r-project.org/web/packages/stabs/stabs.pdf 

## This code is run in parallel on MSI supercomputing nodes.

######################### Initial setup ######################
## before initiating parallel 
rm(list=ls()) ## clear workspace

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel", "data.table") ## package methods is not loaded by default by RScript on MSI
check.packages(packages)

## now initiate parallel processors.
cl <- makeCluster(8) # 8 workers per core -- itasca
registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))

###################### Read the arguments passed from command line in MSI ##########
## For MSI
args <- commandArgs(TRUE)
input.genes.list <- args[1] # "input.dir/genes_split_dir/genes_split_n.txt"
input.genes.table <- args[2] # genes x samples
input.microbiome.table <- args[3] # microbiome x samples
output.dir <- args[4] #output.dir.parallel.allgenes

## debug by spitting output to console
print(paste0("genes.list: ",input.genes.list))
print(paste0("genes.table: ",input.genes.table))
print(paste0("microbiome.table: ",input.microbiome.table))
print(paste0("output.dir: ",output.dir))

######################### Input: Genes and taxa tables ############

genes <- data.frame(fread(input.genes.table,sep="\t",head=T), row.names = 1, check.names = F); dim(genes)
microbes <- data.frame(fread(input.microbiome.table,sep="\t",head=T), row.names = 1, check.names =F); dim(microbes)

print(paste0("# genes = ",dim(genes)[2]))

print(paste0("# microbes = ",dim(microbes)[2]))

## Ensure same sampleIDs in both genes and microbes matrices
stopifnot(all(rownames(genes) == rownames(microbes)))

## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors

##################### Stability selection ####################

stabs_stability <- function( x, y, gene_name){
  
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("stabs","glmnet","methods","doParallel") ## package methods is not loaded by default by RScript. 
  check.packages(packages)
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## perform stability selection using glmnet lasso
  stab.glmnet <- stabsel(x = x, y = y_i,
                         fitfun = glmnet.lasso, cutoff = 0.6,
                         PFER = 1)
  
  taxa.selected <- names(stab.glmnet$selected)
  if(length(taxa.selected) == 0) taxa.selected <-"None"
  
  taxa.selected.df <- as.data.frame(cbind(gene_name,taxa.selected))
  
  return(taxa.selected.df)
  
}

## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")

## invoke parallel computation for fitting model for each gene
parallel_time <- system.time({
  # sink(paste0(output.dir,"/log.txt"), append=TRUE)
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% {
    # cat(paste("Processing gene:",gene_list[i],"\n"))
    # sink() #end diversion of output
    stabs_stability(x,y,gene_list[i])
  }
})
stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_stabs_stabsel_output.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)
