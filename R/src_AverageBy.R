#' Average a numeric matrix or data.frame by columns
#' 
#' Function takes a named list with column numbers to be averaged. Unpaired columns are simply returned as-is.
#' 
#' @param cts a numeric matrix or data.frame
#' @param by a named list with column numbers to be averaged, colname will be the respective list name
#' 
#' @details
#' 
#' @author Alexander Toenges
#' 
#' @examples 
#' cts <- lapply(seq(1,5), function(x) rnorm(10)) %>% do.call(cbind, .)
#' colnames(cts) <- c("A_rep1", "B_rep1", "A_rep2", "C_rep1", "B_rep2")
#' rownames(cts) <- paste0("Gene", seq(1,nrow(cts)))
#' by <- sapply(unique(gsub("_rep.*", "", colnames(cts))), function(x) grep(paste0("^", x, "_rep"), colnames(cts)))
#' AverageBy(cts,by)
#' 
#' @export
AverageBy <- function(cts, by){
  
  invisible(match.arg(class(by), "list"))
  invisible(match.arg(arg = class(cts)[1], choices = c("matrix", "data.frame")))
  if(!length(unlist(by)) == ncol(cts) | max(unlist(by)) > ncol(cts)){
    stop("by has more/less elements than cts has columns", call. = FALSE)
  } 
  if(is.null(names(by))) stop("by must have names", call. = FALSE)
  if(class(cts)[1] != "matrix") cts <- as.matrix(cts)
  
  averaged <- lapply(by, function(x) matrixStats::rowMeans2(cts[,x,drop=FALSE])) %>% 
    do.call(cbind, .)
  
  colnames(averaged) <- names(by)
  rownames(averaged) <- rownames(cts)
  
  return(averaged)
  
}
    
    