## Tutorial
This tutorial demonstrates integration methods used for joint analysis of host transcriptomic and microbiome data as described in _Priya et al. "Shared and disease-specific host gene-microbiome associations across human diseases"_.

### 1. Sparse CCA
Sparse Canonical Correlation Analysis (sparse CCA) identifies linear combination of subsets of variables from two datasets such that they are maximally correlated. We will apply this approach to identify groups of host genes that are associated with groups of microbial taxa.   

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

For further processing to visualize sparse CCA components, perform enrichment analysis on selected genes, etc., please check here (add link).

### 2. Lasso

Demo data for download:
- gene expression data
- microbiome data

Link to script with all functions in Lasso tutorial

Step 1: Read input data
```R
## load gene expression data
genes <- load_gene_expr("gene_expr_demo_lasso.txt")
dim(genes)

## load microbiome data
microbes <- load_microbiome_abnd("microbiome_demo_lasso.txt")
dim(microbes)

## Ensure same samples in both genes and microbes data
stopifnot(all(rownames(genes) == rownames(microbes)))

y <- genes
x <- microbes
```

Step 2: Fit lasso model and perform inference using desparsified lasso
```R
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
```


Step 3: Stability selection

```R
## set a seed for replicability
set.seed(0511)

stab.glmnet <- stabsel(x = x, y = y_i,
                       fitfun = glmnet.lasso, cutoff = 0.6,
                       PFER = 1)

taxa.selected <- names(stab.glmnet$selected)
if(length(taxa.selected) == 0) taxa.selected <-"None"

taxa.selected
```

Step 4: Merge output from 2. and 3. to get associations






