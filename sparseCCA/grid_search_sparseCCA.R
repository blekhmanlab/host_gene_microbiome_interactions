## Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu

## This code implements grid search for sparse CCA analysis


## initial setup 
rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
setwd("/path/to/working_directory")

## import libraries
library(PMA) ##for sparse CCA
library(data.table)

###################### Functions ################
load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F)
  genes <- as.matrix(genes)
  
}

load_microbiome_abnd <- function(filename){
  microbes <- read.table(filename,sep="\t",head=T, row.names = 1, check.names =F)
  microbes <- as.matrix(microbes)
  
}

filter_genes <- function(genes, qt){
  genes.t <- t(genes)
  genes.sd <- transform(as.data.frame(genes.t), SD=apply(as.data.frame(genes.t),1, sd, na.rm = TRUE))
  ## select top genes with high SD (~ variability) across samples
  SD_quantile <- quantile(genes.sd$SD) ## identical to summary(genes.sd$SD)
  SD_cutoff <- SD_quantile[qt] ## 2nd quantile -- 25th quantile.
  genes.sd <- genes.sd[order(genes.sd$SD, decreasing = T),]
  top.variable.genes <- rownames(genes.sd[genes.sd$SD > SD_cutoff,])
  ## subset these genes from gene table
  select <- which(colnames(genes) %in% top.variable.genes)
  genes <- genes[,select]
}

get_tuning_params <- function(X, Z, outputFile = NULL, pdfFile = NULL, num_perm=25){
  perm.out <- CCA.permute(X,Z,typex="standard",typez="standard", nperms = num_perm)
  ## can tweak num. of permutations if needed (default = 25)
  
  ## get p-value for first left and right canonical covariates resulting from the selected
  ## tuning parameter value
  perm.out.pval <- perm.out$pvals[which(perm.out$penaltyxs == perm.out$bestpenaltyx)]
  
  if(!is.null(outputFile)){
    sink(outputFile)
    print(perm.out)
    print(paste0("P-value for selected tuning params:",perm.out.pval)) #0
    sink()
  }
  
  if(!is.null(pdfFile)){
    pdf(pdfFile)
    plot(perm.out)
    dev.off()
  }
  
  return(perm.out)
}

run_sparseCCA <- function(X, Z, CCA.K, penaltyX, penaltyZ, vInit=NULL, outputFile=NULL){
  CCA.out <-  CCA(X,Z,typex="standard",typez="standard",K=CCA.K,
                  penaltyx=penaltyX,penaltyz=penaltyZ,
                  v=vInit) ## standardize=T by default
  if(!is.null(outputFile)){
    sink(outputFile)
    print(CCA.out)
    sink()
  }
  
  ## Inspired by code from Engelhardt's paper: https://github.com/daniel-munro/imageCCA/blob/master/run_CCA.R
  ## add rownames to output factors
  rownames(CCA.out$u) <- colnames(X)
  rownames(CCA.out$v) <- colnames(Z)
  ## compute contribution of selected features to each of the samples.
  CCA_var_genes <- X %*% CCA.out$u ## canonical variance for genes 
  CCA_var_microbes <- Z %*% CCA.out$v ## canonical variance for microbes
  
  return(list(CCA.out, CCA_var_genes, CCA_var_microbes))
  
}

get_avg_features <- function(cca_cov, CCA.K){
  num_features <- 0
  for(k in 1:CCA.K){
    num_features <- num_features + length(which(cca_cov[,k]!=0))
  }
  avg_features <- num_features/CCA.K
}

save_CCA_components <- function(CCA.out, CCA.K, dirname){
  ## Print canonical covariates in files 
  for(i in CCA.K){
    #i <- 2 ##debug
    print(paste0("Writing significant component = ", i))
    selected_X <- which(CCA.out$u[,i]!=0) 
    selected_X <- rownames(CCA.out$u)[selected_X]
    coeff_X <- unname(CCA.out$u[selected_X,i])
    selected_Z <- which(CCA.out$v[,i]!=0)
    selected_Z <- rownames(CCA.out$v)[selected_Z]
    coeff_Z <- unname(CCA.out$v[selected_Z,i])
    ## Make all vectors of same length to avoid repetition of elements from shorter vectors.
    n <- max(length(selected_X), length(selected_Z))
    length(selected_X) <- n                      
    length(selected_Z) <- n
    length(coeff_X) <- n
    length(coeff_Z) <- n
    selected_XZ <- as.data.frame(cbind(gene = selected_X, gene_coeff = coeff_X,
                                       taxa = selected_Z, taxa_coeff = coeff_Z))
    write.table(selected_XZ, file=paste0(dirname,"gene_taxa_component_",i,".txt"), sep = "\t", col.names = NA)
  }
  
}
########### CRC ##############
dataset <- "CRC/case"

## load gene expression and microbiome tables
# genes <- load_gene_expr("./data/clean/CRC/gene_expr/tumor_protein_coding_vsd_t.txt") #[1]    44 16685
genes <- load_gene_expr("./data/clean/CRC/gene_expr/tumor_protein_coding_vsd_SDfilt_t.txt")
dim(genes) #[1]    44 12513
# microbes <- load_microbiome_abnd("./data/clean/CRC/taxa_abnd/tumor_taxa_filt_clr_w_gender.txt") #[1]  44 229 -- old code
# microbes <- load_microbiome_abnd("./data/clean/CRC/taxa_abnd/tumor_taxa_clr_t_no_contam.txt")
microbes <- load_microbiome_abnd("./data/clean/CRC/taxa_abnd/tumor_taxa_clr_t_no_contam_0.001_0.1_V3.txt")
dim(microbes) 
# [1]  44 235

## Ensure same sampleIDs in both genes and microbes data before sparse CCA
stopifnot(all(rownames(genes) == rownames(microbes)))

## select tuning parameters using grid-search
X <- genes
Y <- microbes
scoreXcv <- c()
scoreYcv <- c()
penaltyX <- seq(0.1,0.4,length=10)
penaltyY <- seq(0.15,0.4,length=10)
corr_CRC <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
num_samples <- nrow(genes)
start_time <- Sys.time()
for( i in 1:length(penaltyX)){
  for(j in 1:length(penaltyY)){
    # print(paste0("Index: i = ",i,", j =", j)); flush.console()
    for(k in 1:num_samples){
      ##debug
      # i <- 1
      # j <- 1
      # k <- 2
      ##debug
      print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()
      #compute weights with sample k held out:
      # Default niter = 15 edited to 5 to speed this up.
      res <- CCA(X[-k,],Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F)
      ## Compute scores for k'th sample for first pair of canonical variables
      ## Take weight of features (res$u and res$v) computed using all except 
      ## the kth sample and multiple it by values for the kth sample in the 
      ## feature matrix X and Y. 
      ## %*% implies matrix multiplication. 
      scoreXcv[k] <- X[k,]%*%res$u ## single value
      scoreYcv[k] <- Y[k,]%*%res$v ## single value
    }
    ## correlation between scores for X and Y for all held out samples.
    corr_CRC[i,j] = cor(scoreXcv,scoreYcv) 
  }
}
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste0("Time elapsed for CRC = ", time_elapsed))
# [1] "Time elapsed for CRC = 50.9241799155871"

row.names(corr_CRC) <- as.character(penaltyX)
colnames(corr_CRC) <- as.character(penaltyY)

corr_CRC_df <- as.data.frame(corr_CRC)
rownames(corr_CRC_df)
colnames(corr_CRC_df)

## save to file
# write.table(corr_CRC_df, file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/corr_CRC_0.1_0.4_V4.txt"), sep="\t", row.names = T, col.names = NA )
# save(corr_CRC, file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/corr_CRC_0.1_0.4_V4.RData"))

## load precomputed grid-search output
dataset <- "CRC/case"
penaltyX <- seq(0.1,0.4,length=10)
penaltyY <- seq(0.15,0.4,length=10)
load(paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/corr_CRC_0.1_0.4_V4.RData"))

# find index with max absolute corr
bestpenalty <- which(abs(corr_CRC) == max(abs(corr_CRC)), arr.ind = TRUE)
bestpenalty

bestpenaltyX <- penaltyX[bestpenalty[1]]
bestpenaltyX 
bestpenaltyY <- penaltyY[bestpenalty[2]]
bestpenaltyY

## order abs corr to get top 5 corr
index <- order(abs(corr_CRC), decreasing = T)
abs(corr_CRC)[index][1:5] ## top 5 absolute corr

cca.k = 10 ## number of desired components

## Run sparse CCA using selected tuning param using permutation search
cca <- run_sparseCCA(genes, microbes, cca.k, bestpenaltyX, bestpenaltyY,
                     outputFile=paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/CCA.output.",bestpenaltyX,"_",bestpenaltyY,".txt"))

## canonical correlation for each component:
cca[[1]]$cors

## average number of genes and microbes in resulting components
avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
avg_genes
# [1] 371.6 

avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
avg.microbes
# 12.1

## Test significance of correlation using LOOCV
X <- genes
Y <- microbes
cca.k = 10
scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
corr_pval <- c()
corr_r <- c()
for(i in 1:nrow(genes)){ #n = no. of samples
  #compute weights with sample i held out:
  res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyY, K=cca.k, trace = F) ## default niter = 15 which is spit out when trace = T (default)
  ###compute scores for i'th sample for each component (pair of canonical variables)
  for(j in 1:cca.k){
    print(paste0("i = ", i," K = ", j)); flush.console()
    scoresXcv[i,j] <- X[i,]%*%res$u[,j]
    scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
  }
}
## Test for each components
for(j in 1:cca.k){
  # plot(scoresXcv,scoresYcv)
  corr <- cor.test(scoresXcv[,j],scoresYcv[,j])
  corr_pval[j] <- corr$p.value
  corr_r[j] <- corr$estimate
}
corr_pval

length(which(corr_pval < 0.1)) 
which(corr_pval < 0.1)
# [1]  1  2  3  4  5  6  7  9 10
length(which(corr_pval < 0.05))
which(corr_pval < 0.05)
# [1]  1  2  4  5  6  7  9 10
corr_padj <- p.adjust(corr_pval, method = "BH")
corr_padj
# [1] 3.336897e-03 1.508995e-13 2.169986e-01 1.523028e-02 1.214574e-04 1.923676e-03 2.529803e-09 6.255788e-03 8.735609e-01 6.544138e-01
which(corr_padj < 0.1)
# [1]  1 2 4 5 6 7 8
length(which(corr_padj < 0.1))

## LOOCV corr
corr_r

outputFile <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/gene_taxa_components/crc_sparseCCA_summary_",bestpenaltyX,"_",bestpenaltyY,".txt")
sink(outputFile)
cat(paste0(" bestpenaltyX = ", bestpenaltyX, ", bestpenaltyY = ", bestpenaltyY))
cat(paste0("\n cor(Xu,Yv): \n"))
cat(paste0(signif(cca[[1]]$cors, digits = 4)))
cat(paste0("\n Avg. no. of genes across components = ",avg_genes))
cat(paste0("\n Avg. no. of microbes across components= ", avg.microbes))
cat(paste0("\n P-value for components (LOOCV): \n"))
cat(paste0(signif(corr_pval, digits = 4)))
cat(paste0("\n LOOCV corr: \n"))
cat(paste0(signif(corr_r, digits = 4)))
cat(paste0("\n No. of components with p-value < 0.1 = ", length(which(corr_pval < 0.1))))
cat(paste0("\n No. of components with p-value < 0.05 = ", length(which(corr_pval < 0.05))))
cat(paste0("\n No. of components with FDR < 0.1 = ", length(which(corr_padj < 0.1))))
cat(paste0("\n Significant components: \n" ))
cat(paste0(which(corr_padj < 0.1)))
sink()

# write.table(cca[[2]], file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/CCA_var_genes.txt"), sep="\t", row.names = T, col.names = NA )
# write.table(cca[[3]], file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/CCA_var_microbes.txt"), sep="\t", row.names = T, col.names = NA )

## only spit out significant components
sig <- which(corr_padj < 0.1)
dirname <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/gene_taxa_components/sig_gene_taxa_components_",bestpenaltyX,"_", bestpenaltyY,"_padj/")
## This will return FALSE if the directory already exists or is uncreatable, 
## and TRUE if it didn't exist but was succesfully created.
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
save_CCA_components(cca[[1]],sig,dirname)


########### IBD #############
dataset <- "IBD/case"

## load gene expression and microbiome tables
# genes <- load_gene_expr("./data/clean/IBD/gene_expr/IBD_protein_coding_vsd_t.txt") #[1]    56 15980
genes <- load_gene_expr("./data/clean/IBD/gene_expr/IBD_protein_coding_vsd_SDfilt_t.txt")
dim(genes) #[1]    56 11985
# microbes <- load_microbiome_abnd("./data/clean/IBD/taxa_abnd/IBD_taxa_filt_clr_w_gender.txt") #[1]   56 72
# microbes <- load_microbiome_abnd("./data/clean/IBD/taxa_abnd/IBD_taxa_clr_t.txt")
microbes <- load_microbiome_abnd("./data/clean/IBD/taxa_abnd/IBD_taxa_clr_t_no_contam_0.001_0.1.txt")
dim(microbes) 
#[1]  56 127
# [1]  56 121 -- 3/16/2020, post contam removal. 

## Ensure same sampleIDs in both genes and microbes data before sparse CCA
stopifnot(all(rownames(genes) == rownames(microbes)))

## select tuning parameters using grid-search
X <- genes
Y <- microbes
scoreXcv <- c()
scoreYcv <- c()
penaltyX <- seq(0.1,0.4,length=10)
penaltyY <- seq(0.15,0.4,length=10)
corr_IBD <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
num_samples <- nrow(genes)
start_time <- Sys.time()
for( i in 1:length(penaltyX)){
  for(j in 1:length(penaltyY)){
    # print(paste0("Index: i = ",i,", j =", j)); flush.console()
    for(k in 1:num_samples){
      print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()
      #compute weights with sample k held out:
      # Default niter = 15 edited to 5 to speed this up.
      res <- CCA(X[-k,],Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F)
      ###compute scores for k'th sample for first pair of canonical variables
      scoreXcv[k] <- X[k,]%*%res$u
      scoreYcv[k] <- Y[k,]%*%res$v
    }
    corr_IBD[i,j] = cor(scoreXcv,scoreYcv) ## correlation for all held out samples.
  }
}
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste0("Time elapsed for IBD = ", time_elapsed))
# [1] "Time elapsed for IBD = 1.20009556584888" ## assuming this is in hours.

row.names(corr_IBD) <- as.character(penaltyX)
colnames(corr_IBD) <- as.character(penaltyY)

corr_IBD_df <- as.data.frame(corr_IBD)
rownames(corr_IBD_df)
colnames(corr_IBD_df)
# 
# ## save to file
# write.table(corr_IBD_df, file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/corr_IBD_0.1_0.4_V2.txt"), sep="\t", row.names = T, col.names = NA )
# save(corr_IBD, file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/corr_IBD_0.1_0.4_V2.RData"))
# 

# ## find index with max corr
bestpenalty <- which(abs(corr_IBD_df) == max(abs(corr_IBD_df)), arr.ind = TRUE)
bestpenalty

## load precomputed grid-search output
dataset <- "IBD/case"
penaltyX <- seq(0.1,0.4,length=10)
penaltyY <- seq(0.15,0.4,length=10)
load(paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/corr_IBD_0.1_0.4_V2.RData"))

# ## find index with max corr
bestpenalty <- which(abs(corr_IBD) == max(abs(corr_IBD)), arr.ind = TRUE)
bestpenalty

bestpenaltyX <- penaltyX[bestpenalty[1]]
bestpenaltyX
bestpenaltyY <- penaltyY[bestpenalty[2]]
bestpenaltyY 


cca.k = 10 ## number of desired components

## Run sparse CCA using selected tuning param using permutation search
cca <- run_sparseCCA(genes, microbes, cca.k, bestpenaltyX, bestpenaltyY,
                     outputFile=paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/CCA.output.gridsearch",bestpenaltyX,"_",bestpenaltyY,".txt"))

## canonical correlation for each component:
cca[[1]]$cors

## average number of genes and microbes in resulting components
avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
avg_genes


avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
avg.microbes

## Test significance of correlation using LOOCV 
X <- genes
Y <- microbes
cca.k = 10
scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
corr_pval <- c()
corr_r <- c()
for(i in 1:nrow(genes)){ #n = no. of samples
  #compute weights with sample i held out:
  res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyY, K=cca.k, trace = F) ## default niter = 15 which is spit out when trace = T (default)
  ###compute scores for i'th sample for each component (pair of canonical variables)
  for(j in 1:cca.k){
    print(paste0("i = ", i," K = ", j)); flush.console()
    scoresXcv[i,j] <- X[i,]%*%res$u[,j]
    scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
  }
}
## Test 
for(j in 1:cca.k){
  # plot(scoresXcv,scoresYcv)
  corr <- cor.test(scoresXcv[,j],scoresYcv[,j])
  corr_pval[j] <- corr$p.value
  corr_r[j] <- corr$estimate
}
corr_pval

length(which(corr_pval < 0.1))
which(corr_pval < 0.1)
# [1] 1 2 6 8 

length(which(corr_pval < 0.05))
corr_padj <- p.adjust(corr_pval, method = "BH")
corr_padj
# [1] 5.916186e-04 2.871835e-09 2.782138e-01 4.634485e-01 4.634485e-01 7.641898e-03 2.222044e-01 2.083545e-02 7.129199e-01 2.782138e-01
length(which(corr_padj < 0.1)) 
which(corr_padj < 0.1)
# [1] 1 2 3 4 9

corr_r

outputFile <- paste0("./data/analysis/",dataset,
                     "/output_sparseCCA_V2/grid_search/gene_taxa_components/IBD_case_sparseCCA_",bestpenaltyX,"_",bestpenaltyY,".txt")
sink(outputFile)
cat(paste0(" bestpenaltyX = ", bestpenaltyX, ", bestpenaltyY = ", bestpenaltyY))
cat(paste0("\n cor(Xu,Yv): \n"))
cat(paste0(" ",signif(cca[[1]]$cors, digits = 4)))
cat(paste0("\n Avg. no. of genes across components = ",avg_genes))
cat(paste0("\n Avg. no. of microbes across components= ", avg.microbes))
cat(paste0("\n P-value for components (LOOCV): \n"))
cat(paste0(" ",signif(corr_pval, digits = 4)))
cat(paste0("\n LOOCV corr: \n"))
cat(paste0(" ",signif(corr_r, digits = 4)))
cat(paste0("\n No. of components with p-value < 0.1 = ", length(which(corr_pval < 0.1))))
cat(paste0("\n No. of components with p-value < 0.05 = ", length(which(corr_pval < 0.05))))
cat(paste0("\n No. of components with FDR < 0.1 = ", length(which(corr_padj < 0.1))))
cat(paste0("\n Significant components (FDR p-value < 0.1): \n" ))
cat(paste0(" ",which(corr_padj < 0.1)))
sink()

# write.table(cca[[2]], file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/CCA_var_genes.txt"), sep="\t", row.names = T, col.names = NA )
# write.table(cca[[3]], file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/CCA_var_microbes.txt"), sep="\t", row.names = T, col.names = NA )

## only spit out significant components
sig <- which(corr_padj < 0.1)
dirname <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/gene_taxa_components/sig_gene_taxa_components_",bestpenaltyX,"_", bestpenaltyY,"_padj/")
## This will return FALSE if the directory already exists or is uncreatable, 
## and TRUE if it didn't exist but was succesfully created.
# save_CCA_components(cca[[1]],sig,dirname)


########### IBS ###################

dataset <- "IBS/case"

## load gene expression and microbiome tables
# genes <- load_gene_expr("./data/clean/IBS/gene_expr/IBS_protein_coding_vsd_t.txt") #[1]    29 16568
genes <- load_gene_expr("./data/clean/IBS/gene_expr/IBS_protein_coding_vsd_SDfilt_t.txt")
dim(genes) #[1]    29 12429
# microbes <- load_microbiome_abnd("./data/clean/IBS/taxa_abnd/IBS_taxa_avg_T1_T2_filt_clr_w_gender.txt") #[1]  29 177
# microbes <- load_microbiome_abnd("./data/clean/IBS/taxa_abnd/IBS_taxa_clr_t.txt")
microbes <- load_microbiome_abnd("./data/clean/IBS/taxa_abnd/IBS_taxa_clr_t_no_contam_0.001_0.1.txt")
dim(microbes) 
# [1]  29 238

## Ensure same sampleIDs in both genes and microbes data before sparse CCA
stopifnot(all(rownames(genes) == rownames(microbes)))

## select tuning parameters using grid-search
X <- genes
Y <- microbes
scoreXcv <- c()
scoreYcv <- c()
penaltyX <- seq(0.1,0.4,length=10)
penaltyY <- seq(0.15,0.4,length=10)
corr_IBS <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
num_samples <- nrow(genes)
start_time <- Sys.time()
for( i in 1:length(penaltyX)){
  for(j in 1:length(penaltyY)){
    # print(paste0("Index: i = ",i,", j =", j)); flush.console()
    for(k in 1:num_samples){
      print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()
      #compute weights with sample k held out:
      # Default niter = 15 edited to 5 to speed this up.
      res <- CCA(X[-k,],Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F)
      ###compute scores for k'th sample for first pair of canonical variables
      scoreXcv[k] <- X[k,]%*%res$u
      scoreYcv[k] <- Y[k,]%*%res$v
    }
    corr_IBS[i,j] = cor(scoreXcv,scoreYcv) ## correlation for all held out samples.
  }
}
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste0("Time elapsed for IBS = ", time_elapsed))
# [1] "Time elapsed for IBS = 26.7878766496976"

row.names(corr_IBS) <- as.character(penaltyX)
colnames(corr_IBS) <- as.character(penaltyY)

corr_IBS_df <- as.data.frame(corr_IBS)
rownames(corr_IBS_df)
colnames(corr_IBS_df)
 
# ## save to file
# write.table(corr_IBS_df, file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/corr_IBS_0.1_0.4.txt"), sep="\t", row.names = T, col.names = NA )
# save(corr_IBS, file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/corr_IBS_0.1_0.4.RData"))

## load precomputed grid-search output
# dataset <- "IBS/case"
# penaltyX <- seq(0.1,0.4,length=10)
# penaltyY <- seq(0.15,0.4,length=10)
# load(paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/corr_IBS_0.1_0.4.RData"))

# ## find index with max corr
bestpenalty <- which(corr_IBS == max(corr_IBS), arr.ind = TRUE)
bestpenalty

## find index with max absolute corr
bestpenalty <- which(abs(corr_IBS) == max(abs(corr_IBS)), arr.ind = TRUE)
bestpenalty

## Go with max abs error, high LOOCV bsolute value -0.8713614
bestpenaltyX <- penaltyX[bestpenalty[1]] 
bestpenaltyX
bestpenaltyY <- penaltyY[bestpenalty[2]]
bestpenaltyY 

cca.k = 10 ## number of desired components

## Run sparse CCA using selected tuning param using permutation search
cca <- run_sparseCCA(genes, microbes, cca.k, bestpenaltyX, bestpenaltyY,
                     outputFile=paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/CCA.output.gridsearch_",bestpenaltyX,"_",bestpenaltyY,".txt"))

## canonical correlation for each component:
cca[[1]]$cors

## average number of genes and microbes in resulting components
avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
avg_genes


avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
avg.microbes

## Test significance of correlation using LOOCV
X <- genes
Y <- microbes
cca.k = 10
scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
corr_pval <- c()
corr_r <- c()
for(i in 1:nrow(genes)){ #n = no. of samples
  #compute weights with sample i held out:
  res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyY, K=cca.k, trace = F) ## default niter = 15 which is spit out when trace = T (default)
  ###compute scores for i'th sample for each component (pair of canonical variables)
  for(j in 1:cca.k){
    print(paste0("i = ", i," K = ", j)); flush.console()
    scoresXcv[i,j] <- X[i,]%*%res$u[,j]
    scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
  }
}
## Test 
for(j in 1:cca.k){
  # plot(scoresXcv,scoresYcv)
  corr <- cor.test(scoresXcv[,j],scoresYcv[,j])
  corr_pval[j] <- corr$p.value
  corr_r[j] <- corr$estimate
}
corr_pval

## Note, not significant for some components after 5th component (pval < 0.1) 
length(which(corr_pval < 0.1)) 
which(corr_pval < 0.1)
# [1]  1  2  3  5  9 10

length(which(corr_pval < 0.05)) 
corr_padj <- p.adjust(corr_pval, method = "BH")
corr_padj
# [1] 8.052846e-10 2.289964e-04 2.592038e-04 8.628416e-01 3.881318e-03 2.780899e-01 8.628416e-01 3.111814e-01 2.735869e-02 5.138806e-02
length(which(corr_padj < 0.1)) 
which(corr_padj < 0.1)


corr_r

outputFile <- paste0("./data/analysis/",dataset,
                     "/output_sparseCCA_V2/grid_search/gene_taxa_components/IBS_case_sparseCCA_",bestpenaltyX,"_",bestpenaltyY,".txt")
sink(outputFile)
cat(paste0(" bestpenaltyX = ", bestpenaltyX, ", bestpenaltyY = ", bestpenaltyY))
cat(paste0("\n cor(Xu,Yv): \n"))
cat(paste0(" ",signif(cca[[1]]$cors, digits = 4)))
cat(paste0("\n Avg. no. of genes across components = ",avg_genes))
cat(paste0("\n Avg. no. of microbes across components= ", avg.microbes))
cat(paste0("\n P-value for components (LOOCV): \n"))
cat(paste0(" ",signif(corr_pval, digits = 4)))
cat(paste0("\n LOOCV corr: \n"))
cat(paste0(" ",signif(corr_r, digits = 4)))
cat(paste0("\n No. of components with p-value < 0.1 = ", length(which(corr_pval < 0.1))))
cat(paste0("\n No. of components with p-value < 0.05 = ", length(which(corr_pval < 0.05))))
cat(paste0("\n No. of components with FDR < 0.1 = ", length(which(corr_padj < 0.1))))
cat(paste0("\n Significant components (FDR p-value < 0.1): \n" ))
cat(paste0(" ",which(corr_padj < 0.1)))
sink()

# write.table(cca[[2]], file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/CCA_var_genes.txt"), sep="\t", row.names = T, col.names = NA )
# write.table(cca[[3]], file = paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/CCA_var_microbes.txt"), sep="\t", row.names = T, col.names = NA )
sig <- which(corr_padj < 0.1)
dirname <- paste0("./data/analysis/",dataset,"/output_sparseCCA_V2/grid_search/gene_taxa_components/sig_gene_taxa_components_",bestpenaltyX,"_", bestpenaltyY,"_padj/")
## This will return FALSE if the directory already exists or is uncreatable, 
## and TRUE if it didn't exist but was succesfully created.
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
save_CCA_components(cca[[1]],sig,dirname)


