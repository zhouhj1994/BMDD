# BMDD
Bimodal Dirichlet distribution (BMDD) for microbiome data

Reference: Huijuan Zhou, Jun Chen, and Xianyang Zhang. (2023). Bimodal Dirichlet Distributions And Its Application to Microbiome Data Analysis.

The package implements a variational EM algorithm for estimating the parameters in bimodal Dirichlet distribution.
For the microbiome data, the unobserved true compositions are assumed to follow the bimodal Dirichlet distribution,
thus the package also estimates the posterior distribution of the true compositions.

## Installation
```r
devtools::install_github("zhouhj1994/BMDD")
```
