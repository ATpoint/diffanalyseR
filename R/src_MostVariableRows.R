#' Get most variable rows from a numeric matrix#' 
#' 
#' Find most variable rows of a numeric matrix. If rownames are present returns those,
#' if not returns row index.
#' 
#' @param counts numeric matrix
#' @param ntop number of most most variable rows to return
#' @param type c("rv", "mgv") rv: selection by rowVars and mgv: selection based on 
#' the scran function \code{modelGeneVar} followed by \code{getTopHVGs} (take the top ntop) and return rownames
#' 
#' @details 
#'  
#' @author Alexander Toenges
#' 
#' @examples
#' counts <-lapply(seq(1,4), function(x) rnorm(20)) %>% do.call(cbind, .)
#' MostVariableRows(counts, 10)
#' 
#' @export
MostVariableRows <- function(counts, ntop = 1000, type = c("rv", "mgv")){

  invisible(match.arg(arg = class(counts)[1], choices = c("matrix")))
  invisible(match.arg(arg = class(ntop), choices = c("numeric", "integer")))
  type <- match.arg(type)
  
  if(type=="rv"){
    ## ntop most variable genes/regions:    
    rv <- matrixStats::rowVars(as.matrix(counts))
    
    if(length(rv) < ntop) message("Fewer than ntop genes in counts. Returning all genes.")
    
    selected <- head(order(rv, decreasing = TRUE), ntop)
  
    if(is.null(rownames(counts))){
      
      message("No rownames found in matrix, returning index of top variable rows")
      return(selected)
      
    } else return(rownames(counts)[selected])
  }
  
  if(type=="mgv"){
    
    if(is.null(rownames(counts))) {
      rwn=FALSE; message("No rownames found in matrix, returning index of top variable rows")
    }
    
    return(scran::getTopHVGs(scran::modelGeneVar(x=counts),n=ntop, row.names=rwn))
    
  }
  
}
