## Author: Sambhawa Priya
## Blekhman Lab
## priya030@umn.edu 

## This code implements gene-wise lasso model to identify set of gut microbes 
## associated with host gene expression. We use desparsified lasso for inference
## and 95% CI for host gene-microbe associations. 

rm(list=ls()) ## ## check to make sure no precomputed dataframe before clearing

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel") ## package methods is not loaded by default by RScript. 
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

## Spit output to console
print(paste0("genes.list: ",input.genes.list))
print(paste0("genes.table: ",input.genes.table))
print(paste0("microbiome.table: ",input.microbiome.table))
print(paste0("output.dir: ",output.dir))

################# Input genes and taxa matrix ###########

## could be  generalized further by passing as arg. For now good for debugging.
# data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, comment.char = "")
genes <- read.table(input.genes.table,sep="\t",head=T, row.names = 1, check.names = F); dim(genes)
microbes <- read.table(input.microbiome.table,sep="\t",head=T, row.names = 1, check.names =F); dim(microbes)

## Transform genes to make samples as rows and features as columns
# genes <- t(genes)
# microbes <- t(microbes) ## these already have samples as rows

## Ensure same sampleIDs in both genes and microbes matrices
stopifnot(all(rownames(genes) == rownames(microbes)))

## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors

## In downstream analysis, estimateSigma() throws error with duplicated columns in x. 
## Hence let's fix it here
# any(duplicated(x, MARGIN = 2))
# # [1] TRUE
# which(duplicated(x, MARGIN = 2))
# 45 
## Figure out later why this is duplicated.
## [1] "Bacteria;Bacteroidetes;Flavobacteriia;Flavobacteriales"

x.uniq <- unique(x, MARGIN = 2)
dim(x.uniq)
# [1]  44 228
x <- x.uniq

print(paste0("# genes = ",dim(y)[2]))

print(paste0("# microbes = ",dim(x)[2]))

############### Functions ############

estimate.sigma <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    
    sigmas <- numeric()
    ## Repeat 5x and take median of finite values else assign sigma to ss_res and raise the flag 
    for( j in 1:5){
      bestlambda = get.lambda(x, y_i, 10, 10)
      lasso.fit = glmnet(x,y_i,alpha = 1)
      y = as.vector(y_i)
      yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
      beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
      df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
      n = length(y_i)
      ss_res = sum((y - yhat)^2)
      if((n-df-1) >= 1) {
        sigmas[j] = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
      }
    }
    if(length(sigmas != 0 )){
      
      sigma = median(sigmas, na.rm = T)
      sigma.flag = 1
    } else{
      sigma = 1 ## conservative option
      # sigma = ss_res ## lenient option
      sigma.flag = 2
    }
    
  }
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}

## Edit: Aug 10, 2019: Inspired from the details for estimateSigma() function in selectiveInference package: https://cran.r-project.org/web/packages/selectiveInference/selectiveInference.pdf 
## Seems very close to the computation of sigma in hdi lasso package: do.initial.fit() in https://github.com/cran/hdi/blob/master/R/helpers.R
estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
      sigma = 1 ## conservative option
      # sigma = ss_res ## lenient option
      sigma.flag = 2
    }
    
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}


fit.cv.lasso <- function(x, y_i, kfold, repeats){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  
  for (i in 1:repeats){
    
      ## glmnet CV
      fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE)  
      errors = data.frame(fit$lambda,fit$cvm)
      lambdas <- rbind(lambdas,errors)
      
      ## Get R^2 of final model
      r.sqr.final[i] <- r_squared(as.vector(y_i), 
                                  as.vector(predict(fit$glmnet.fit, 
                                                    newx = x, s = fit$lambda.min)))
      ## Get adjusted R^2
      r.sqr.final.adj[i] <- adj_r_squared(r.sqr.final[i], n = nrow(x), 
                                          p = sum(as.vector(coef(fit$glmnet.fit, 
                                                                 s = fit$lambda.min)) > 0))
      ## cross-validated test-set R^2
      ## Adapted from: https://rpubs.com/kaz_yos/alasso
      # test_R2 <- lapply(unique(fit$foldid), function(id) {
      #   ## Fit excluding test set (foldid == id)
      #   # print(id); flush.console()
      #   fit <- glmnet(x = x[fit$foldid != id,],
      #                 y = y_i[fit$foldid != id],
      #                 alpha = 1,
      #                 lambda = fit$lambda.min) ## glmnet manual advises against proving single lambda value.  
      #   ## Test-set Y_hat using model fit at best lambda
      #   y_pred <- predict(fit, newx = x[fit$foldid == id,], s = fit$lambda.min)
      #   ## Test-set Y
      #   y <- y_i[fit$foldid == id]
      #   ## Test-set R^2
      #   r_sqr <- 1 - sum((y - y_pred)^2) / sum((y - mean(y))^2)
      #   # print(r_sqr)
      #   return(r_sqr)
      # })
      # 
      # r.sqr.CV.test[i] <- median(unlist(test_R2))
  }
  
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)
  # dim(lambdas)
  # select the best one
  bestindex = which(lambdas[2]==min(lambdas[2]))
  bestlambda = lambdas[bestindex,1]
  
  return(list(bestlambda = bestlambda, r.sqr = median(r.sqr.final), 
              r.sqr.adj = median(r.sqr.final.adj) #, r.sqr.CV.test = median(r.sqr.CV.test)
              ))
}

## functions to compute R2
## Borrowed from: https://rpubs.com/kaz_yos/alasso
r_squared <- function(y, yhat) {
  ybar <- mean(y)
  ## Total SS
  ss_tot <- sum((y - ybar)^2)
  ## Residual SS
  ss_res <- sum((y - yhat)^2)
  ## R^2 = 1 - ss_res/ ss_tot
  1 - (ss_res / ss_tot)
}
## Function for Adjusted R^2
## n sample size, p number of prameters
adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}

############### Estimate sigma (standard deviation of the error term or noise) ###############

estimate.sigma.fit.hdi <- function(x, y, gene_name){
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("glmnet","hdi","methods","doParallel") ## package methods is not loaded by default by RScript. 
  check.packages(packages)
  
  ## Debug for local test
  # library(glmnet)
  # library(hdi)
  # gene_name <- "SNX21" ## debug for a specific gene
  
  print(paste0("Processing gene:", gene_name));flush.console()
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## Fit lasso CV model
  fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i), repeats = 1)
  bestlambda <- fit.model$bestlambda
  r.sqr <- fit.model$r.sqr
  r.sqr.adj <- fit.model$r.sqr.adj
  # r.sqr.CV.test <- fit.model$r.sqr.CV.test
 
  ## Estimate sigma using the estimated lambda param
  sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
  sigma <- sigma.myfun$sigmahat
  beta <- as.vector(sigma.myfun$betahat)[-1]
  sigma.flag <- sigma.myfun$sigmaflag
  
  ## Edit: Aug 10. Add Filter here to compute inference for only the genes with R2 > 0.25. This is our current pipeline. 
  
  ## Inference 
  lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
  lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
  lasso.FDR.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)), 
                             taxa = names(lasso.proj.fit$pval.corr), 
                             r.sqr = r.sqr, r.sqr.adj = r.sqr.adj,
                             pval = lasso.proj.fit$pval, padj = lasso.proj.fit$pval.corr, 
                             ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                             sigma = sigma, sigma.flag = sigma.flag,
                             row.names=NULL)
  
  ##local debug
  # View(lasso.FDR.df)
  
  ## Edit Oct 23, 2019 -- No need to save the intermediate RData files. 
  # save(lasso.FDR.df, file=paste0(output.dir,"/",gene_name,"_lasso_hdi.RData")) 
  
  return(lasso.FDR.df)
  
}

## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")

## invoke parallel computation for fitting model for each gene
parallel_time <- system.time({
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% estimate.sigma.fit.hdi(x,y,gene_list[i])
})
stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

## Print result for all the genes processed 
# print(parallel_res)

## Print result for this node's genes list
filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_lasso_hdi.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)

