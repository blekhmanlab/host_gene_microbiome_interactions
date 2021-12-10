# Host gene-microbiome associations
A framework in R for identifying associations between gut microbiome composition and host gene expression in human diseases. These scripts are used in the study Priya et al. "Shared and disease-specific host gene-microbiome associations across human diseases".

### Description of scripts in sparse CCA pipeline:

#### grid_search_sparseCCA.R
This script performs Sparse Canonical Correlation Analysis (sparse CCA) to identify groups of host genes that are associated with groups of microbial taxa in each disease cohort. R package PMA is used to compute sparse CCA components. The script first performs hyperparameter tuning using grid-search, then fits the sparse CCA model using the tuned paramters, and outputs significant sparse CCA components. 

#### msigdb_enrichment_sparseCCA.R
This script performs enrichment analysis for host genes found associated with gut microbes in sparse CCA components. Erichment analysis using MSigDB pathway database is performed per component, followed by correction for multiple hypothesis testing (FDR). Pathways are pooled across components, sorted by FDR q-value, and duplicated pathways are filtered. This also includes comparitive log odds-ratio based differential enrichment analysis to compare host pathways between case and control groups within a disease cohort.

#### collapse_redundant_pathways_sparseCCA.R
This script computes similarity metrics between pathways to identify and collapse redundant pathways for visualization.

#### overlap_enriched_pathways_sparseCCA.R
This script identifies overlaps between host pathways across disease cohorts to identify disease-specific and common pathways that are associated with gut microbiome across diseases.

#### overlap_components_across_diseases_sparseCCA.R
This script computes overlapping genes between sparse CCA components across different disease cohorts. It also performs enrichment analysis to identify host pathways enriched among genes shared across sparse CCA components that associate with gut microbes across diseases.

