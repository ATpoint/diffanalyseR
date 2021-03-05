#' Column-wise binning of a numeric matrix
#' 
#' Perform binning on a numeric matrix. 
#' Divide the matrix in ncol(mat)/byNcols columns and then on each of these
#' apply rowMeans, combine into a binned output matrix.
#' 
#' @param mat numeric matrix
#' @param byNcols average this number of consecutive columns
#' @param method aggregation method, "mean" or "median"
#' @param drop.uneven logical, whether to drop columns if ncol(mat)/byNcols is odd, see details.
#' 
#' @details 
#' `drop.uneven` determines what to do with columns if ncol(mat)/byNcols is odd.
#' Say you have 99 columns and want to bin into 10er columns, means the last 9 would be orphans.
#' With FALSE these would simply be averaged and returned with the rest, or TRUE would be discarded.
#' If mat has fewer columns that one sets byNcols then with FALSE the function returns all columns averaged,
#' or with TRUE a stop is returned.
#' 
#' @author Alexander Toenges
#' 
#' @examples 
#' mat <- sapply(paste0("Sample",seq(1,99)), function(x) rnorm(100,20,3))
#' binned <- binMatrix(mat = mat, byNcols = 3, method = "mean")
#' @export
binMatrix <- function(mat, byNcols, method = "mean", drop.uneven = FALSE){
  
  invisible(match.arg(arg = class(mat)[1], choices = c("matrix")))
  invisible(match.arg(arg = class(byNcols)[1], choices = c("numeric")))
  invisible(match.arg(arg = method, choices = c("mean", "median")))
  
  if(byNcols==1) return(mat)
  
  full <- floor(ncol(mat) / byNcols)
  half <- ncol(mat) - byNcols*full
  
  if(byNcols > ncol(mat)){
    if(drop.uneven) {
      stop("byNcols > ncol(mat) and drop.uneven=FALSE", call. = FALSE)
    } else {
      warning("byNcols > ncol(mat) aggregating all columns into a single one", call. = FALSE)
      full<-1; half<-0; byNcols<-ncol(mat)
    }
  }
  
  #/ the columns to bin over:
  op <- lapply(seq(1, full*byNcols, byNcols), function(x) seq(x,x+byNcols-1))
  if(half!=0 & drop.uneven == "FALSE") {
    k <- max(op[[length(op)]])
    op[[length(op)+1]] <- seq(k+1, k+half)
  }
  
  names(op) <- paste0("bin", seq(1, length(op)))
  
  if(method=="mean") func=matrixStats::rowMeans2
  if(method=="median") func=matrixStats::rowMedians
  
  return(sapply(op, function(x) func(mat[,x, drop=FALSE])))
  
}
  




  
  
  
  
