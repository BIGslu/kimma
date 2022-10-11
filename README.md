# kimma

[![DOI](https://zenodo.org/badge/387951897.svg)](https://zenodo.org/badge/latestdoi/387951897)

Kinship In Mixed Model Analysis

Dill-McFarland KA, Mitchell K, Batchu S, Segnitz RM, Benson B, Janczyk T, Cox MS, Mayanja-Kizza HA, Boom WH, Benchek P, Stein C, Hawn TR, Altman MC. 2022. kimma: Flexible linear mixed effects modeling with kinship for RNA-seq data. bioRxiv doi: [10.1101/2022.10.10.508946](https://doi.org/10.1101/2022.10.10.508946)

We introduce kimma (Kinship In Mixed Model Analysis), an open-source R package for flexible linear mixed effects modeling of RNA-seq including covariates, weights, random effects, covariance matrices, and fit metrics. In simulated datasets, kimma detects differentially expressed genes (DEGs) with similar specificity, sensitivity, and computational time as limma unpaired and dream paired models. Unlike other software, kimma supports covariance matrices as well as fit metrics like AIC. Utilizing genetic kinship covariance, kimma revealed that kinship impacts model fit and DEG detection in a related cohort. Thus, kimma equals or outcompetes current DEG pipelines in sensitivity, computational time, and model complexity. 

### Installation

```
install.packages("devtools")
devtools::install_github("BIGslu/kimma")
```

# Vignette

<https://bigslu.github.io/kimma_vignette/kimma_vignette.html>

Other related [vignettes](https://bigslu.github.io/tutorials/) and [workshops](https://bigslu.github.io/workshops/)
