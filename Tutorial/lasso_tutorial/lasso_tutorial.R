## Lasso tutorial
## Sambhawa Priya

## Goal
## Our aim in this tutorial is to show how to run our lasso analysis on a small set of genes (~2-3 genes)
## to identify host gene-taxa associations for the demo dataset.


## Install and import libraries
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("glmnet","data.table","hdi","stabs") 
check.packages(packages)

############## Functions ###############

## Please execute these functions

load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
  genes <- as.matrix(genes)
  
}

load_microbiome_abnd <- function(filename){
  microbes <- read.table(filename,sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
  microbes <- as.matrix(microbes)
  
}

estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  

  ## Fit a lasso object
  lasso.fit = glmnet(x,y_i,alpha = 1) ## this is same as cv.fit$glmnet.fit from loocv code below.
  beta <- as.vector(coef(lasso.fit, s = bestlambda)) ## This gives coefficients of fitted model, not predicted coeff.
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  ## predicted coefficients, same as coefficient of fitted model lasso. Either one is fine.
  # beta = predict(lasso.fit,s=bestlambda, type="coef")
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1))
    sigma.flag = 0
  } else{
    sigma = 1 ## conservative option
    # sigma = ss_res ## lenient option
    sigma.flag = 2
  }
  
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}

fit.cv.lasso <- function(x, y_i, kfold){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  
  ## glmnet CV
  cv.fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE, standardize = T)  
  lambdas = data.frame(cv.fit$lambda,cv.fit$cvm)
  
  ## get best lambda -- lambda that gives min cvm
  bestlambda <- cv.fit$lambda.min
  bestlambda_index <- which(cv.fit$lambda == bestlambda)
  
  ## Get R^2 of final model
  final_model <- cv.fit$glmnet.fit
  r_sqr_final_model <- cv.fit$glmnet.fit$dev.ratio[bestlambda_index]
  
  ## Get adjusted R^2
  r_sqr_final_adj <- adj_r_squared(r_sqr_final_model, n = nrow(x), 
                                   p = sum(as.vector(coef(cv.fit$glmnet.fit, 
                                                          s = cv.fit$lambda.min)) > 0))
  
  return(list(bestlambda = bestlambda, r.sqr = r_sqr_final_model, 
              r.sqr.adj = r_sqr_final_adj
  ))
}

## functions to compute R2
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


############### Input lasso demo data #############

# ## In Rstudio, find the path to the directory where the current script is located.
# current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# ## current_dir should point to path of "lasso_tutorial".
# 
# 
# genes <- load_gene_expr(paste0(current_dir,"/input/gene_expr_demo_lasso.txt"))
# dim(genes) #[1]    44 3
# 
# ## load microbiome data (note, here we load microbiome data with sex covariate)
# microbes <- load_microbiome_abnd(paste0(current_dir,"/input/microbiome_demo_lasso.txt"))
# dim(microbes) #[1]  44 236 -- 235 taxa + 1 sex covariate
# 
# 
# ## Ensure same sampleIDs in both genes and microbes data before sparse CCA
# stopifnot(all(rownames(genes) == rownames(microbes)))
# 
# y <- genes #response
# x <- microbes #predictors
# 
# 
# ############ Fit lasso model and test inference using HDI ############
# 
# ## We are going to test three genes: WNT5A, RIPK3, and SMAP2 for their association with microbes
# 
# ## Extract expression of first gene in the matrix
# i <- 3 ## replace with 2 or 3 to test other two genes
# y_i <- y[,i]
# gene_name <- colnames(y)[i]
# 
# ## Make sure y_i is numeric before model fitting
# stopifnot(class(y_i) == "numeric")
# 
# ## Fit lasso CV model
# fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i))
# bestlambda <- fit.model$bestlambda
# r.sqr <- fit.model$r.sqr ## note this will give us R^2 for the gene's final model fit using bestLambda
# ## This R^2 reflects final model R^2 for this gene using all the microbes in the model,
# ## and does not correspond to each gene-microbe pair.
# 
# ## Estimate sigma and betainit using the estimated LOOCV lambda.
# ## Sigma is the standard deviation of the error term or noise.
# sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
# sigma <- sigma.myfun$sigmahat
# beta <- as.vector(sigma.myfun$betahat)[-1] ## remove intercept term
# sigma.flag <- sigma.myfun$sigmaflag
# 
# ## Inference using lasso projection method, also known as the de-sparsified Lasso,
# ## using an asymptotic gaussian approximation to the distribution of the estimator.
# lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
# ## A few lines of log messages appear here along with a warning about substituting sigma value (standard deviation of error term or noise)
# ## because we substituted value of sigma using our computation above.
# # Warning message:
# #   Overriding the error variance estimate with your own value.
# 
# ## get 95% confidence interval (CI)
# lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
# 
# ## prep lasso output dataframe
# lasso.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)),
#                        taxa = names(lasso.proj.fit$pval.corr),
#                        r.sqr = r.sqr,
#                        pval = lasso.proj.fit$pval,
#                        ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
#                        row.names=NULL)
# 
# 
# ## sort by p-value
# lasso.df <- lasso.df[order(lasso.df$pval),]
# head(lasso.df)
# 
# ################# Stability selection #################
# 
# ## set a seed for replicability
# set.seed(0511)
# 
# ## perform stability selection using glmnet lasso
# stab.glmnet <- stabsel(x = x, y = y_i,
#                        fitfun = glmnet.lasso, cutoff = 0.6,
#                        PFER = 1)
# 
# taxa.selected <- names(stab.glmnet$selected)
# if(length(taxa.selected) == 0) taxa.selected <-"None"
# 
# 
# stabsel.df <- data.frame("gene" = gene_name, "taxa" = taxa.selected)
# if(taxa.selected == "none"){
#   stabsel.df$stability_selected = "no"
# }else stabsel.df$stability_selected = "yes"
# 
# head(stabsel.df)
# 
# ################ Merge output of lasso+hdi and stabsel #################
# 
# overlap_lasso_stabsel <- merge(lasso.df,stabsel.df, by = c("gene","taxa"))
# head(overlap_lasso_stabsel)

## Explanation of the output
# For first two genes at index 1 and 2 in y (i.e. WNT5A and RIPK3), a taxa is stability selected,
## however, for the 3rd gene, no taxa is stability selected, hence we have an empty dataframe after merging
## outputs of lasso and stability selection.

## In step 2. "Fit lasso model and test inference using desparsified lasso", you can toggle index for
## y between 1, 2, and 3 to test the pipeline for different genes.
