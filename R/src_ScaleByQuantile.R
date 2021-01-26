#' Winsorize a count matrix using fixed quantiles by row
#' 
#' Per row trim the values beyond the specified quantiles to that quantile,
#' can be useful e.g. when plotting heatmaps to limit influence of outliers
#' Takes a numeric matrix or data.frame and returns one with the winsorized values.
#' 
#' @param Counts a matrix or dataframe with the count data
#' @param lower bottom quantile
#' @param upper top quantile
#' @author Alexander Toenges
#' 
#' @examples 
#' 
#' Counts <- sapply(seq(1,10), function(x) rnorm(1000,100))
#' ScaleByQuantile(Counts=Counts,lower=0.1,upper=0.9)
#'  
#' @export
ScaleByQuantile <- function(Counts, lower = 0, upper = 1){
  
  if(class(Counts)[1] != "matrix" & class(Counts)[1] != "data.frame"){
    stop("Counts must be a matrix or data.frame")
  }
  
  if(lower==0 & upper==1) return(Counts)
  
  if(class(Counts)[1] == "data.frame"){
    cnames <- colnames(Counts)
    Counts <- as.matrix(counts)
    colnames(Counts) <- cnames
  }
  
  if(upper < 1){
    qt.upper <- as.numeric(quantile(Counts, upper, na.rm = TRUE))
    Counts[Counts > qt.upper] <- qt.upper
  }
  if(lower > 0){
    qt.lower <- as.numeric(quantile(Counts, lower, na.rm = TRUE))
    Counts[Counts < qt.lower] <- qt.lower
  }
  
  return(Counts)
  
}
