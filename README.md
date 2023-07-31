# BMDD
Bimodal Dirichlet distribution (BMDD) for microbiome data

Reference: Huijuan Zhou, Lu Yang, Jun Chen, and Xianyang Zhang. (2023). Bimodal Dirichlet Distributions And Its Application to Microbiome Data Analysis.

The package implements a variational EM algorithm for estimating the parameters in bimodal Dirichlet distribution.
For the microbiome data, the unobserved true compositions are assumed to follow the bimodal Dirichlet distribution,
thus the package also estimates the posterior distribution of the true compositions.

## Installation
```r
devtools::install_github("zhouhj1994/BMDD")
```
## Example
This is an example of using BMDD as an zero-imputation method for LinDA to identify colorectal cancer associated bacterial species. Data "phy" is a phyloseq-class experiment-level object. LinDA is a differential abundance analysis method for microbiome data. The package is available at https://CRAN.R-project.org/package=MicrobiomeStat and https://github.com/zhouhj1994/LinDA.

Reference for the dataset: Yu et al. (2017). [Metagenomic analysis of faecal microbiome as a tool towards targeted non-invasive biomarkers for colorectal cancer](https://gut.bmj.com/content/66/1/70)
Reference for LinDA: Zhou et al. (2022). [LinDA: linear models for differential abundance analysis of microbiome compositional data] (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02655-5)

```r
#install package "phyloseq" for importing "phy" dataset
data(phy)
otu_filter <- function(feature.dat, prev = 0.1, dep = 1000){
  idx <- apply(feature.dat, 1, function(x) sum(x > 0) > (ncol(feature.dat) * prev))
  idx2 <- colSums(feature.dat) > dep
  return(feature.dat[idx, idx2])
}
otu.tab <- as.data.frame(as.matrix(otu_table(phy)))
meta.dat <- as.data.frame(as.matrix(sample_data(phy)))
meta.dat$grp <- as.factor(meta.dat$grp)
feature.dat <- otu_filter(otu.tab)
meta.dat <- meta.dat[colnames(feature.dat), ]

bmdd.fit <- bmdd(W = feature.dat, type = 'count')
prop.bmdd <- t(t(bmdd.fit$beta) / colSums(bmdd.fit$beta))
bmdd.obj  <- MicrobiomeStat::linda(feature.dat = prop.bmdd, meta.dat = meta.dat,
                                   formula = '~grp', feature.dat.type = 'proportion')
bmdd.res <- bmdd.obj$output[[1]][,'padj',drop = F]

linda.obj  <- MicrobiomeStat::linda(feature.dat = feature.dat, meta.dat = meta.dat,
                                    formula = '~grp', feature.dat.type = 'count')
linda.res <- linda.obj$output[[1]][,'padj',drop = F]
```
