## Tutorial

### 1. Sparse CCA

Demo data for download: 
- gene expression data
- microbiome data

Link to script with all functions in sparse CCA tutorial

Step 1: Read input data

```R
## load gene expression data
genes <- load_gene_expr("gene_expresion_demo.txt")
dim(genes)

## load microbiome data
microbes <- load_microbiome_abnd("microbiome_demo.txt")
dim(microbes)

## Ensure same samples in both genes and microbes data
stopifnot(all(rownames(genes) == rownames(microbes)))
```

Step 2: Tune hyperparameters

This step uses grid-search that takes a while, so we precomputed penalty values for demo dataset that you can set and skip to Step 3.
```R
bestpenaltyX <- 0.05
bestpenaltyY <- 0.3222
```

```R
## Skip this if using pre-computed values
## select tuning parameters using grid-search
bestPenalty <- tune_params_grid_search(genes,microbes)
bestpenaltyX <- bestPenalty[1]
bestpenaltyY <- bestPenalty[2]
```

Step 3: Run sparse CCA

```R
## Set the number of desired CCA components
cca.k = 10

## Create a folder named "sparseCCA_output_demo" in your current working directory for this demo.
cca <- run_sparseCCA(genes, microbes, cca.k, bestpenaltyX, bestpenaltyY,
                     outputFile=paste0("./sparseCCA_output_demo/CCA.output.",bestpenaltyX,"_",bestpenaltyY,".txt"))

## average number of genes and microbes in resulting components
avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
avg_genes

avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
avg.microbes
```

Step 4: Test significance of components using LOOCV

```R
CCA_pval <- test_significance_LOOCV(genes, microbes, bestpenaltyX, bestpenaltyY, cca.k)

length(which(CCA_pval < 0.1)) 
which(CCA_pval < 0.1)

CCA_padj <- p.adjust(CCA_pval, method = "BH")
CCA_padj

length(which(CCA_padj < 0.1))
which(CCA_padj < 0.1)
```

Step 5. Output significant components at FDR < 0.1

```R
sig <- which(CCA_padj < 0.1)
dirname <- paste0("./sparseCCA_output_demo/gene_taxa_components/sig_gene_taxa_components_",bestpenaltyX,"_", bestpenaltyY,"_padj/")
## This returns False if the directory already exists or can't be created, 
## returns True if it didn't exist but was successfully created.
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
save_CCA_components(cca[[1]],sig,dirname)
```

### 2. Lasso

Demo data for download:
- gene expression data
- microbiome data

Link to script with all functions in Lasso tutorial

Step 1: Read input data

Step 2: Fit lasso model and perform inference using HDI

Step 3: Stability selection

Step 4: Merge output from 2. and 3. to get associations






