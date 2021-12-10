# Host gene-microbiome associations
A framework in R for identifying associations between gut microbiome composition and host gene expression in human diseases. These scripts are used in the study Priya et al. "Shared and disease-specific host gene-microbiome associations across human diseases".

### Description of scripts used in the sparse CCA pipeline:

#### - grid_search_sparseCCA.R
This script performs Sparse Canonical Correlation Analysis (sparse CCA) to identify groups of host genes that are associated with groups of microbial taxa in each disease cohort. R package _PMA_ is used to compute sparse CCA components. The script first performs hyperparameter tuning using grid-search, then fits the sparse CCA model using the tuned paramters, and outputs significant sparse CCA components. 

#### - msigdb_enrichment_sparseCCA.R
This script performs enrichment analysis for host genes found associated with gut microbes in sparse CCA components. Erichment analysis using MSigDB pathway database is performed per component, followed by correction for multiple hypothesis testing (FDR). Pathways are pooled across components, sorted by FDR q-value, and duplicated pathways are filtered. This also includes comparitive log odds-ratio based differential enrichment analysis to compare host pathways between case and control groups within a disease cohort.

#### - collapse_redundant_pathways_sparseCCA.R
This script computes similarity metrics between pathways to identify and collapse redundant pathways for visualization.

#### - overlap_enriched_pathways_sparseCCA.R
This script identifies overlaps between host pathways across disease cohorts to identify disease-specific and common pathways that are associated with gut microbiome across diseases.

#### - overlap_components_across_diseases_sparseCCA.R
This script computes overlapping genes between sparse CCA components across different disease cohorts. It also performs enrichment analysis to identify host pathways enriched among genes shared across sparse CCA components that associate with gut microbes across diseases.

### Description of scripts used in the lasso pipeline:

#### - hdi_lasso_proj_loocv_MSI.R
This script implements a lasso regression for each host gene’s expression as response and abundances of microbial taxa and values of other covariates as independent variables. It uses leave-one-out cross-validation to estimate the tuning parameter, which is used to fit the final model on a given disease dataset. It uses desparsified lasso approach (R package _HDI_) to obtain 95% confidence intervals and p-values for the coefficient of each microbe associated with a given host gene. This script implements a parallel framework for executing the gene-wise lasso analysis, where execution of lasso models on host genes can be parallelized across multiple nodes and cores on a compute cluster.

#### - stabs_stability_selection.R
This script performs stability selection (using R package _stabs_) for the lasso model to identify microbes robustly associated with host genes. This script also uses parallel execution framework described above.

#### - postprocess_lasso_stabsel_output_V2.R
This script performs postprocessing of lasso output and stability selection runs to identify significant and stability selected host gene-taxa associations.

#### - identify_case_specific_associations_features.R
This script identifies case-specific associations and features across disease cohorts by factoring out control associations/features from case associations/features.

#### - enrichment_lasso_disease_specific_genes.R
This script performs enrichment analysis for host genes found significantly associated with a taxa in each disease.

#### - overlap_CRC_IBD_IBS.R
This script identifies overlaps between host gene-taxa associations across disease cohorts. It identifies gut microbial taxa and host genes that are shared between associations across diseases, and identifies their network of associations in each disease cohort for visualization.
