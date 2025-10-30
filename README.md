# BMDD
Bimodal Dirichlet distribution (BMDD) for microbiome data

Reference: Huijuan Zhou, Jun Chen, and Xianyang Zhang. (2025).[BMDD: A probabilistic framework for accurate imputation of zero-inflated microbiome sequencing data](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013124)

The package implements a variational EM algorithm for estimating the parameters in bimodal Dirichlet distribution.
For the microbiome data, the unobserved true compositions are assumed to follow the bimodal Dirichlet distribution.
The package also estimates the posterior distribution of the true compositions.

## Installation
```r
devtools::install_github("zhouhj1994/BMDD")
```
## Example
This is an example of using BMDD as an zero-imputation method for LinDA to identify colorectal cancer associated bacterial species. Data "phy" is a phyloseq-class experiment-level object. LinDA is a differential abundance analysis method for microbiome data. The package is available at https://CRAN.R-project.org/package=MicrobiomeStat and https://github.com/zhouhj1994/LinDA.

- Reference for the dataset: Yu et al. (2017). [Metagenomic analysis of faecal microbiome as a tool towards targeted non-invasive biomarkers for colorectal cancer](https://gut.bmj.com/content/66/1/70)
- Reference for LinDA: Zhou et al. (2022). [LinDA: linear models for differential abundance analysis of microbiome compositional data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02655-5)

```r
#install packages "phyloseq", "MicrobiomeStat"
library(BMDD)

data(phy)

otu_filter <- function(feature.dat, prev = 0.1, dep = 1000){
  idx <- apply(feature.dat, 1, function(x) sum(x > 0) > (ncol(feature.dat) * prev))
  idx2 <- colSums(feature.dat) > dep
  return(feature.dat[idx, idx2])
}
feature.dat <- as.data.frame(as.matrix(phyloseq::otu_table(phy)))
meta.dat <- as.data.frame(as.matrix(phyloseq::sample_data(phy)))
meta.dat$grp <- as.factor(meta.dat$grp)
feature.dat <- otu_filter(feature.dat)
meta.dat <- meta.dat[colnames(feature.dat), ]

m <- nrow(feature.dat)
n <- ncol(feature.dat)

bmdd.obj <- bmdd(W = feature.dat, type = 'count', trace = TRUE)

# posterior mean
beta <- bmdd.obj$beta
post.mean <- t(t(beta) / colSums(beta))

## generate 100 posterior samples of the composition for each sample
zero.fun <- function(X) {
  X <- t(apply(X, 1, function (x) {
    if(all(x == 0)) {
      x[x == 0] <- min(X[X != 0])
    } else {
      x[x == 0] <- min(x[x != 0]) 
    }
    return(x)
  }))
  return(X)
}

K <- 100
beta <- beta[, rep(1 : n, K)]
X <- matrix(rgamma(m * n * K, beta, 1), m)
X <- t(t(X) / colSums(X))
if(any(X == 0)) {
  X <- zero.fun(X)
  X <- t(t(X) / colSums(X))
}
colnames(X) <- paste0('sample', 1 : (n * K))
rownames(X) <- rownames(beta)

## apply LinDA to the proportion matrix, with 100 replicates per sample.
id <- factor(rep(1 : n, K))
grp <- rep(meta.dat$grp, K)
Z <- cbind.data.frame(grp, id)
rownames(Z) <- colnames(X)

linda.bmdd.obj <- MicrobiomeStat::linda(feature.dat = X, meta.dat = Z, 
                                   formula = "~grp+(1|id)", feature.dat.type = "proportion")

# apply LinDA to the original count matrix
linda.obj <- MicrobiomeStat::linda(feature.dat = feature.dat, meta.dat = meta.dat, 
                                   formula = "~grp", feature.dat.type = "count")

```
