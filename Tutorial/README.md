## Tutorial
This tutorial demonstrates integration methods used for joint analysis of host transcriptomic and microbiome data as described in _Priya et al. "Shared and disease-specific host gene-microbiome associations across human diseases"_.

### 1. Sparse CCA
Sparse Canonical Correlation Analysis (sparse CCA) identifies linear combination of subsets of variables from two datasets such that they are maximally correlated (Witten et al. 2009). We will apply this approach to identify groups of host genes that are associated with groups of microbial taxa.   

#### Setting up input and output

Download the folder [_sparseCCA_tutorial_](https://github.com/blekhmanlab/host_gene_microbiome_interactions/blob/main/Tutorial/sparseCCA_tutorial.zip) at a relevant location on your computer. This folder includes the Rscript _sparseCCA_tutorial.R_ with all functions used in sparse CCA tutorial, an input folder with demo dataset (i.e. _gene_expresion_demo_sp_CCA.txt_ and _microbiome_demo_sp_CCA.txt_), and an output folder. 
Open and execute the script _sparseCCA_tutorial.R_ to load all libraries and functions, and follow the steps below.    

#### Step 1: Read input data

```R
## In Rstudio, find the path to the directory where the current script is located.
## If not using Rstudio, current_dir should point to your working directory for this demo.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

## load gene expression data
genes <- load_gene_expr(paste0(current_dir,"/input/gene_expresion_demo_sp_CCA.txt"))
dim(genes)

## load microbiome data
microbes <- load_microbiome_abnd(paste0(current_dir,"/input/microbiome_demo_sp_CCA.txt"))
dim(microbes)

## Ensure same samples in both genes and microbes data
stopifnot(all(rownames(genes) == rownames(microbes)))
```

#### Step 2: Tune hyperparameters

This step uses grid-search that takes a while to run, so we've pre-computed penalty values for the demo dataset that you can set as follows, and skip to Step 3.
```R
bestpenaltyX <- 0.05
bestpenaltyY <- 0.3222
```

```R
## Skip to Step 3 if using pre-computed values above
## select tuning parameters using grid-search
bestPenalty <- tune_params_grid_search(genes,microbes)
bestpenaltyX <- bestPenalty[1]
bestpenaltyY <- bestPenalty[2]
```

#### Step 3: Run sparse CCA

```R
## Set the number of desired CCA components
cca.k = 10

## run sparse CCA
cca <- run_sparseCCA(genes, microbes, cca.k, bestpenaltyX, bestpenaltyY,
                     outputFile=paste0(current_dir,"/output/CCA_demo_output_",bestpenaltyX,"_",bestpenaltyY,".txt"))

## average number of genes and microbes in resulting components
avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
avg_genes

avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
avg.microbes
```

#### Step 4: Test significance of components using leave-one-out cross-validation

```R
## This will take ~1 min to run. 
CCA_pval <- test_significance_LOOCV(genes, microbes, bestpenaltyX, bestpenaltyY, cca.k)

## which components have p-value < 0.1
length(which(CCA_pval < 0.1)) 
which(CCA_pval < 0.1)

## adjust for multiple testing
CCA_padj <- p.adjust(CCA_pval, method = "BH")
length(which(CCA_padj < 0.1))
which(CCA_padj < 0.1)
```

#### Step 5. Output significant components

```R
sig_cutoff <- 0.1 
sig <- which(CCA_padj < sig_cutoff)
dirname <- dirname <- paste0(current_dir,"/output/demo_gene_taxa_components/")
## This returns returns true if directory didn't exist but was successfully created,
## and returns false if the directory already exists or can't be created.
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
save_CCA_components(cca[[1]],sig,dirname)
```
Each output sparse CCA component includes a set of host genes that are correlated with a set of taxa. It also includes non-zero weights (or canonical loadings) on gut microbes, and non-zero weights on a subset of host genes correlated with those gut microbes to capture joint variation in the two sets of observations.   

For further processing to visualize sparse CCA components, perform enrichment analysis on selected genes, etc., please check [here](https://github.com/blekhmanlab/host_gene_microbiome_interactions/tree/main/sparseCCA).

### 2. Lasso

Lasso is a penalized regression approach that uses shrinkage or regularization to perform variable selection (Tibshirani et al. 1996). We will use lasso penalized regression to identify specific associations between individual host genes and gut microbial taxa. To identify microbial taxa that are correlated with a host gene, we model host gene as response and abundances of microbiome taxa as independent variables. 

#### Setting up input and output

Download the folder [_lasso_tutorial.zip_](https://github.com/blekhmanlab/host_gene_microbiome_interactions/blob/main/Tutorial/lasso_tutorial.zip) at a relevant location on your computer. This folder includes the Rscript _lasso_tutorial.R_ with all functions used in the lasso tutorial and an input folder with demo dataset (i.e. _gene_expresion_demo_lasso.txt_ and _microbiome_demo_lasso.txt_). 
Open and execute the script _lasso_tutorial.R_ to load all libraries and functions, and follow the steps below.

#### Step 1: Read input data
```R
## In Rstudio, find the path to the directory where the current script is located.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

## load gene expression data
genes <- load_gene_expr(paste0(current_dir,"/input/gene_expr_demo_lasso.txt"))
dim(genes)

## load microbiome data
microbes <- load_microbiome_abnd(paste0(current_dir,"/input/microbiome_demo_lasso.txt"))
dim(microbes)

## Ensure same samples in both genes and microbes data
stopifnot(all(rownames(genes) == rownames(microbes)))

y <- genes #response
x <- microbes #predictors
```

#### Step 2: Fit lasso model and perform inference using desparsified lasso
```R
## We are going to test three genes: WNT5A, RIPK3, and SMAP2 for their association with microbes

## Extract expression of first gene in the matrix
i <- 1 ## replace with 2 or 3 to test other two genes
y_i <- y[,i]
gene_name <- colnames(y)[i]

## Make sure y_i is numeric before model fitting
stopifnot(class(y_i) == "numeric")

## Fit lasso model using LOOCV
fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i))
bestlambda <- fit.model$bestlambda
r.sqr <- fit.model$r.sqr

## Estimate sigma and betainit using the estimated LOOCV lambda.
## Sigma is the standard deviation of the error term or noise.
sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
sigma <- sigma.myfun$sigmahat
beta <- as.vector(sigma.myfun$betahat)[-1] ## remove intercept term

## Perform inference using lasso projection method, also known as the de-sparsified Lasso,
## using an asymptotic gaussian approximation to the distribution of the estimator.
lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
## A few lines of log messages appear here along with a warning about substituting sigma value (i.e. standard deviation 
## of error term or noise) because we substituted value of sigma using our computation above.
## Warning message:
##   Overriding the error variance estimate with your own value.


## get 95% confidence interval for gene-taxa association
lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))

## prep lasso output dataframe
lasso.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)),
                       taxa = names(lasso.proj.fit$pval.corr),
                       r.sqr = r.sqr,
                       pval = lasso.proj.fit$pval,
                       ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                       row.names=NULL)

## sort by p-value
lasso.df <- lasso.df[order(lasso.df$pval),]
head(lasso.df)
```

#### Step 3: Stability selection

```R
## set a seed for replicability
set.seed(0511)

## perform stability selection using glmnet lasso
stab.glmnet <- stabsel(x = x, y = y_i,
                       fitfun = glmnet.lasso, cutoff = 0.6,
                       PFER = 1)

taxa.selected <- names(stab.glmnet$selected)
if(length(taxa.selected) == 0) taxa.selected <-"None"

stabsel.df <- data.frame("gene" = gene_name, "taxa" = taxa.selected)
if(taxa.selected == "none"){
  stabsel.df$stability_selected = "no"
}else stabsel.df$stability_selected = "yes"

head(stabsel.df)
```

#### Step 4: Merge output from steps 2. and 3. to get gene-taxa associations
```R
overlap_lasso_stabsel <- merge(lasso.df,stabsel.df, by = c("gene","taxa"))
head(overlap_lasso_stabsel)
```

In step 2. "Fit lasso model and test inference using desparsified lasso", you can toggle index for
y between 1, 2, and 3 to test the pipeline for different genes.

For first two genes at index 1 and 2 in y (i.e. _WNT5A_ and _RIPK3_), a taxa is stability selected, however, for the third gene at index 3 (_SMAP2_), no taxa is stability selected, hence we have an empty dataframe after merging outputs of lasso and stability selection.


For further processing of lasso output, please check [here](https://github.com/blekhmanlab/host_gene_microbiome_interactions/tree/main/lasso).

Please note that the parameters used in this tutorial are for the demo dataset used here. These parameters might need to be tweaked for different datasets.

