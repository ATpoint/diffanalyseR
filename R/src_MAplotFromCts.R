#' Pairwise MA-plots using smoothScatter for QC purposes
#' 
#' Given two vectors of normalized counts make an MA-plot to see how normalization performed
#' and whether things look as expected.
#' 
#' @param x counts of sample1
#' @param y counts of sample2
#' @param cts.arelog2 logical, whether counts are already on log2-scale
#' @param type the type of the plot, currently only smoothScatter from base is available, later a ggplot version will come
#' @param abline0 logical, whether to add a horizontal ablone at y=0
#' @param abline2 numeric, add horizontal ablines at +/- abline2 value
#' @param use.loess logical, whether to add a y ~ x  loess fitted line
#' @param ... further arguments passed to smoothMAplot()
#' 
#' @author Alexander Toenges
#' 
#' @details the plot will be directly plotted to devise, so best would be to capture function directly to disk
#' 
#' @examples
#' dds  <- estimateSizeFactors(makeExampleDESeqDataSet())
#' x <- counts(dds,normalized=TRUE)[,1]
#' y <- counts(dds,normalized=TRUE)[,2]
#' MAplotFromCts(x = x,y = y, use.loess=TRUE, Y.limits=c(-2,2))
#' @export
MAplotFromCts <- function(x, y, 
                          Plot.title = "",
                          cts.arelog2 = FALSE,
                          type = "base",
                          abline0 = FALSE,
                          abline2 = NULL,
                          use.loess = FALSE,
                          ...){
  
  #invisible(match.arg(arg = type, choices = c("base")))
  #invisible(match.arg(arg = as.character(cts.arelog2), choices = c("TRUE", "FALSE")))
  #invisible(match.arg(arg = abline0, choices = c("TRUE", "FALSE")))
  #invisible(match.arg(arg = class(abline2), choices = c("numeric")))
  #invisible(match.arg(arg = use.loess, choices = c("TRUE", "FALSE")))
  
  #--------------------------------------------------------------------------------
  
  #/ We always add a prior of 1 to avoid division by zero:
  x <- x + 1
  y <- y + 1
  
  if(cts.arelog2){
    
    logFC    <- x - y
    baseMean <- 0.5 * (x + y)
    
  } else {
    
    logFC    <- log2(x / y)
    baseMean <- 0.5 * log2(x * y)
    
  }
  
  #--------------------------------------------------------------------------------
  
  #/ Plot base::smoothScatter:
  if(type == "base"){
    
    smoothMAplot(X.value = baseMean, 
                 Y.value = logFC, 
                 use.loessFit = use.loess,
                 abline0 = abline0, abline2 = abline2, 
                 Plot.title = Plot.title, 
                 Plot.title.size = 1.5,
                 ...)
    
  }
  
}
  
