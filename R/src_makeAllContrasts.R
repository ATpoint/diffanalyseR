#' Make all possible unique pairwise contrasts
#' 
#' Expecting a vector with group names, e.g. `c("group1", "group2", "group3")` this function will
#' then create all unique pairwise contrasts. 
#' By default the output is a matrix/array in the limma-ish style of `makeContrasts()`.
#' If `deseq2=TRUE` it returns a named list with DESeq2-ish contrasts that can be used with the contrasts argument
#' of `results()` and `lfcShrink()`.
#' The group vector should not contain hyphen!
#' 
#' @param group vector with group names
#' @param delim a string to be used as delimiter for the names, Default is "_vs_", e.g. group1_vs_group2
#' @param deseq2 logical, whether to return contrasts as a name list in DESeq2-ish style, like c("group", "A", "B")
#' @param name if deseq2, then the "group" name
#' 
#' @details 
#' 
#' @author Alexander Toenges
#' 
#' @examples 
#' makeAllContrasts(group = c("A", "B", "C"))
#' makeAllContrasts(group = c("A", "B", "C"), deseq2 = TRUE)
#' 
#' @export
makeAllContrasts <- function(group, 
                             delim = "_vs_", 
                             deseq2 = FALSE,
                             name = "group"){
  
  group <- sort(unique(as.character(group)))
  if(sum(grepl("-", unlist(group)))>0) stop("There is a hyphen somewhere in group", call. = FALSE)
  
  cb <- combn(group, 2)
  
  Contrasts <- list()
  for(x in seq(1, ncol(cb))) Contrasts[[x]] <- paste0(cb[1,x],"-",cb[2,x])
  Contrasts <- limma::makeContrasts(contrasts = unlist(Contrasts), levels = group)
  colnames(Contrasts) <- gsub("-",delim, colnames(Contrasts)) 
  
  message("Created ", ncol(Contrasts), " contrasts")
  
  if(deseq2){
    
    return(sapply(colnames(Contrasts), function(x){
      
      sp <- strsplit(x, split = delim)[[1]]
      return(c(name, sp[1], sp[2]))
      
    }, simplify = FALSE, USE.NAMES = TRUE))
    
  } else return(Contrasts)
  
}
