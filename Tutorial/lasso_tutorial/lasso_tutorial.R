## Lasso tutorial script

## import libraries
library(glmnet)
library(hdi)
library(stabs)

############## Functions ###############

load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
  genes <- as.matrix(genes)
  
}

load_microbiome_abnd <- function(filename){
  microbes <- read.table(filename,sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
  microbes <- as.matrix(microbes)
  
}

adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}

fit.cv.lasso <- function(x, y_i, kfold){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  
  ## glmnet CV
  cv.fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE)  
  lambdas = data.frame(cv.fit$lambda,cv.fit$cvm)
  
  ## get best lambda -- lambda that gives min cvm
  bestlambda <- cv.fit$lambda.min
  bestlambda_index <- which(cv.fit$lambda == bestlambda)
  
  ## Get R^2 of full model (fitted at lambda.min)
  final_model <- cv.fit$glmnet.fit
  r_sqr_final_model <- cv.fit$glmnet.fit$dev.ratio[bestlambda_index]
  # r.sqr.final <- r_squared(as.vector(y_i), 
  #                             as.vector(predict(fit$glmnet.fit, 
  #                                               newx = x, s = fit$lambda.min)))
  # all.equal(r.sqr.final, r_sqr_final_model) ## all.equal uses some toleranceå
  # # [1] TRUE
  
  ## Get adjusted R^2
  r_sqr_final_adj <- adj_r_squared(r_sqr_final_model, n = nrow(x), 
                                   p = sum(as.vector(coef(fit$glmnet.fit, 
                                                          s = fit$lambda.min)) > 0))
  
  return(list(bestlambda = bestlambda, r.sqr = r_sqr_final_model, 
              r.sqr.adj = r_sqr_final_adj))
}


estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  
  ## Our implementation matches do.initial.fit() in hdi code:
  ## https://github.com/cran/hdi/blob/master/R/helpers.R.
  ## The only difference is they compute optimal lambda using 10-fold CV, 
  ## while we compute lambda using LOOCV that we have passed to this function as "bestlambda"
  ## Our logic for this function also matches description for estimateSigma() from selectiveInference package:
  ## https://cran.r-project.org/web/packages/selectiveInference/selectiveInference.pdf
  ## "A lasso regression is fit, using cross-validation to estimate the tuning parameter lambda. With sample size
  ## n, yhat equal to the predicted values and df being the number of nonzero coefficients from the
  ## lasso fit, the estimate of sigma is sqrt(sum((y-yhat)^2) / (n-df-1))."
  
  ## Fit a lasso object
  lasso.fit = glmnet(x,y_i,alpha = 1) ## this is same as cv.fit$glmnet.fit from loocv code below.
  beta <- as.vector(coef(lasso.fit, s = bestlambda)) 
  
  y = as.vector(y_i)
  ## same as hdi code used for computing residual.vector
  ## residual.vector <- y-predict(glmnetfit,newx=x,s=lambda).
  ## In hdi, lambda was computed using 10 fold cv
  ## while here we use lambda computed from our LOOCV code (bestlambda).
  
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1))
  } else{
    sigma = 1 ## conservative option
  }
  
  
  return(list(sigmahat = sigma, betahat = beta)) ## we return beta to be used later in hdi function.
  
}


############# Run lasso to get associations ##########

#### 1: Read data

## load gene expression data
genes <- load_gene_expr("gene_expr_demo.txt")
dim(genes)

## load microbiome data
microbes <- load_microbiome_abnd("microbiome_demo.txt")
dim(microbes)

## Ensure same samples in both genes and microbes data
stopifnot(all(rownames(genes) == rownames(microbes)))

y <- genes
x <- microbes

#### 2: Fit lasso model using LOOCV and perform inference

## Extract the expression for each gene
y_i <- y[,i]

## Make sure y_i is numeric before model fitting 
stopifnot(class(y_i) == "numeric")

## Fit lasso CV model
fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i))
bestlambda <- fit.model$bestlambda
r.sqr <- fit.model$r.sqr

## Estimate sigma using the estimated lambda param
sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
sigma <- sigma.myfun$sigmahat
beta <- as.vector(sigma.myfun$betahat)[-1]

## Inference 
lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
## Throws warning because we substituted our computed sigma (standard deviation of error term or noise)
# Warning message:
#   Overriding the error variance estimate with your own value.

## get 95% CI
lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))

## collect all metrics
lasso.FDR.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)), 
                           taxa = names(lasso.proj.fit$pval), 
                           full_model_r_sqr = r.sqr,
                           pval = lasso.proj.fit$pval, 
                           ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                           row.names=NULL)

#### 3: Stability selection
## set a seed for replicability
set.seed(0511)

## perform stability selection using glmnet lasso
stab.glmnet <- stabsel(x = x, y = y_i,
                       fitfun = glmnet.lasso, cutoff = 0.6,
                       PFER = 1)

taxa.selected <- names(stab.glmnet$selected)
if(length(taxa.selected) == 0) taxa.selected <-"None"

taxa.selected

#### 4. Overlap results from Step 2 and 3.







