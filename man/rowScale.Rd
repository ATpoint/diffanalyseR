% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/src_rowScale.R
\name{rowScale}
\alias{rowScale}
\title{Row scaling with matrixStats}
\usage{
rowScale(
  x,
  center = TRUE,
  scale = TRUE,
  add_attr = FALSE,
  rows = NULL,
  cols = NULL
)
}
\arguments{
\item{x}{input matrix}

\item{center}{logical, whether to center}

\item{scale}{logical, whether to scale}

\item{add_attr}{logical, whether to add scale/center the attributes}

\item{rows}{rows to subset for}

\item{cols}{cols to subset for}
}
\description{
An efficient function for row scaling using matrixStats, see reference for origin.
}
\details{

}
\examples{

}
\references{
A Faster Scale Function (2016) https://www.r-bloggers.com/a-faster-scale-function/
}
\author{
Alexander Toenges
}
