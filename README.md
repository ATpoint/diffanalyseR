# diffanalyseR

Function to automate common tasks related to analysis of NGS count data in R.

## Installation

```r

#/ Check for BiocManager:
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install()
}

#/ Bioc packages:
pk <- c("edgeR", "GenomicRanges")
BiocManager::install(pk)

install.packages("remotes")
remotes::install_github("ATpoint/diffanalyseR")

```

## Brief overview over the functions

See the help of each functions for details and examples.  

- `AverageBy` averages a numeric matrix by the indices of a named list  

- `binMatrix` bins a matrix, for example a 20-col matrix with 4 bins would mean that 5 consecutive columns get averaged into one  

- `gosty` is a wrapper around `gprofiler2::gost` for functional enrichment of gene lists  

- `makeAllContrasts` takes a character of comparison groups and then forms all unique pairwise contrasts in either edgeR or DESeq2 style  

- `MAplotFromCts` produces MAplots with `base::smoothScatter` based on a normalized count matrix, mainly to explore whether sample have been properly normalized  

- `MostVariableRows` finds the most variable rows based on rowwise variance  

- `NormOffsetsFromTxLength` calculate log-offsets that can be used in edgeR and company to correct for both sequencing depth and composition (e.g. via TMM) plus corrects for average transcript length (e.g. via tximport)  

- `Permut_GRanges` calculate permutation-based p-values comparing an observed overlap between two GRanges objects with the expected overlap between one of the objects and a background set  

- `res_edgeR` automates contrast testing via the edgeR QLF framework  

- `rowScale` is a matrixStats-based wrapper to efficiently scale and center numeric matrices as efficient replacement of `t(scale(t(x)))` avoiding the double transposition which takes time with very large data such as single-cell stuff   

- `ScaleByQuantile` winsorizes a numeric matrix based on user-defined quantiles. This can be useful when plotting heatmaps to avoid skewed color breaks due to outliers  