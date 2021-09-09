## Sparse CCA tutorial
## Sambhawa Priya

## Install and import libraries
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("PMA","data.table") 
check.packages(packages)

################# Functions ################

load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
  genes <- as.matrix(genes)
  
}

load_microbiome_abnd <- function(filename){
  microbes <- read.table(filename,sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
  microbes <- as.matrix(microbes)
  
}

run_sparseCCA <- function(X, Z, CCA.K, penaltyX, penaltyZ, vInit=NULL, outputFile=NULL){
  CCA.out <-  CCA(X,Z,typex="standard",typez="standard",K=CCA.K,
                  penaltyx=penaltyX,penaltyz=penaltyZ,
                  v=vInit)
  if(!is.null(outputFile)){
    sink(outputFile)
    print(CCA.out)
    sink()
  }
  
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

tune_params_grid_search <- function( X, Y){
  
  scoreXcv <- c()
  scoreYcv <- c()
  penaltyX <- seq(0.05,0.4,length=10)
  penaltyY <- seq(0.05,0.4,length=10)
  corr_demo <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
  num_samples <- nrow(genes)
  start_time <- Sys.time()
  for( i in 1:length(penaltyX)){
    for(j in 1:length(penaltyY)){
      # print(paste0("Index: i = ",i,", j =", j)); flush.console()
      for(k in 1:num_samples){
       
        print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()
        #compute weights with sample k held out:
        #Default niter = 15 edited to 5 to speed this up.
        res <- CCA(X[-k,],Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F, standardize = T)
        ## Compute scores for k'th sample for first pair of canonical variables
        ## Take weight of features (res$u and res$v) computed using all except 
        ## the kth sample and multiple it by values for the kth sample in the 
        ## feature matrix X and Y.
        scoreXcv[k] <- X[k,]%*%res$u ## single value
        scoreYcv[k] <- Y[k,]%*%res$v ## single value
      }
      ## correlation between scores for X and Y for all held out samples.
      corr_demo[i,j] = cor(scoreXcv,scoreYcv) 
    }
  }
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste0("Time elapsed for param tuning = ", time_elapsed)); flush.console()
  
  row.names(corr_demo) <- as.character(penaltyX)
  colnames(corr_demo) <- as.character(penaltyY)
  
  corr_demo_df <- as.data.frame(corr_demo)
  rownames(corr_demo_df)
  colnames(corr_demo_df)
  
  ##identify best penalty parameters
  # find index with max absolute corr
  bestpenalty <- which(abs(corr_demo) == max(abs(corr_demo)), arr.ind = TRUE)
  bestpenalty
  bestpenaltyX <- penaltyX[bestpenalty[1]]
  bestpenaltyY <- penaltyY[bestpenalty[2]]
  
  return (c(bestpenaltyX,bestpenaltyY))
}

test_significance_LOOCV <- function(X, Y, bestpenaltyX, bestpenaltyY, num_components){
  cca.k = num_components
  scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
  scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
  corr_pval <- c()
  corr_r <- c()
  for(i in 1:nrow(genes)){ #n = no. of samples
    #compute weights with sample i held out:
    res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyY, K=cca.k, trace = F) ## default niter = 15 which is spit out when trace = T (default)
    ###compute scores for i'th sample for each component (pair of canonical variables)
    for(j in 1:cca.k){
      #print(paste0("i = ", i," K = ", j)); flush.console()
      scoresXcv[i,j] <- X[i,]%*%res$u[,j]
      scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
    }
  }
  ## Test for each components
  for(j in 1:cca.k){
    corr <- cor.test(scoresXcv[,j],scoresYcv[,j]) ## Pearson correlation.
    corr_pval[j] <- corr$p.value
  }
  corr_pval
}

########## Tune and run sparse CCA #########

# ## In Rstudio, find the path to the directory where the current script is located.
# current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# 
# #### load data
# 
# ## load gene expression data
# genes <- load_gene_expr(paste0(current_dir,"/input/gene_expresion_demo_sp_CCA.txt"))
# dim(genes)
# 
# ## load microbiome data
# microbes <- load_microbiome_abnd(paste0(current_dir,"/input/microbiome_demo_sp_CCA.txt"))
# dim(microbes)
# 
# ## Ensure same sampleIDs in both genes and microbes data before sparse CCA
# stopifnot(all(rownames(genes) == rownames(microbes)))
# 
# ## set penalty parameters
# bestpenaltyX <- 0.05
# bestpenaltyY <- 0.3222
# 
# ## SKIP if using pre-computed values above
# ## select tuning parameters
# # bestPenalty <- tune_params_grid_search(genes,microbes)
# # bestpenaltyX <- bestPenalty[1]
# # bestpenaltyY <- bestPenalty[2]
# 
# #### Run sparse CCA
# 
# ## Set the number of desired components
# cca.k = 10
# 
# ## Run sparse CCA using selected tuning param using permutation search
# cca <- run_sparseCCA(genes, microbes, cca.k, bestpenaltyX, bestpenaltyY,
#                      outputFile=paste0(current_dir,"/output/CCA_demo_output_",bestpenaltyX,"_",bestpenaltyY,".txt"))
# 
# ## average number of genes and microbes in resulting components
# avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
# avg_genes
# 
# avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
# avg.microbes
# 
# #### Test significance of components using LOOCV
# CCA_pval <- test_significance_LOOCV(genes, microbes, bestpenaltyX, bestpenaltyY, cca.k)
# 
# length(which(CCA_pval < 0.1))
# which(CCA_pval < 0.1)
# 
# CCA_padj <- p.adjust(CCA_pval, method = "BH")
# CCA_padj
# 
# length(which(CCA_padj < 0.1))
# which(CCA_padj < 0.1)
# 
# #### Output significant components
# sig_cutoff <- 0.1
# sig <- which(CCA_padj < sig_cutoff)
# dirname <- paste0(current_dir,"/output/demo_gene_taxa_components/")
# ## This will return FALSE if the directory already exists or is uncreatable,
# ## and TRUE if it didn't exist but was succesfully created.
# ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
# save_CCA_components(cca[[1]],sig,dirname)




