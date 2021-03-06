% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/src_binMatrix.R
\name{binMatrix}
\alias{binMatrix}
\title{Column-wise binning of a numeric matrix}
\usage{
binMatrix(mat, byNcols, method = "mean", drop.uneven = FALSE)
}
\arguments{
\item{mat}{numeric matrix}

\item{byNcols}{average this number of consecutive columns}

\item{method}{aggregation method, "mean" or "median"}

\item{drop.uneven}{logical, whether to drop columns if ncol(mat)/byNcols is odd, see details.}
}
\description{
Perform binning on a numeric matrix. 
Divide the matrix in ncol(mat)/byNcols columns and then on each of these
apply rowMeans, combine into a binned output matrix.
}
\details{
`drop.uneven` determines what to do with columns if ncol(mat)/byNcols is odd.
Say you have 99 columns and want to bin into 10er columns, means the last 9 would be orphans.
With FALSE these would simply be averaged and returned with the rest, or TRUE would be discarded.
If mat has fewer columns that one sets byNcols then with FALSE the function returns all columns averaged,
or with TRUE a stop is returned.
}
\examples{
mat <- sapply(paste0("Sample",seq(1,99)), function(x) rnorm(100,20,3))
binned <- binMatrix(mat = mat, byNcols = 3, method = "mean")
}
\author{
Alexander Toenges
}
