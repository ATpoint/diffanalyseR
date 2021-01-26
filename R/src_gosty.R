#' Functional enrichment of a gene list with gprofiler2
#' 
#' Function minimally expects a vector of genes which will be used for functional enrichment analysis
#' with GOST (g:profiler2). Optionally (and recommended) is a list of background genes.
#' One can specify organism and databases to scan against. See `?gost`, this here is just a simple wrapper to bring
#' output into a more parsable format.
#' 
#' @param InputGenes vector with gene names (HGNC/MGI/Ensembl)
#' @param BackgroundGenes vector with gene names (HGNC/MGI/Ensembl)
#' @param Sources the databases to use, default is c("KEGG", "REAC")
#' @param Species mmusculus/hsapiens etc...
#' @param FDR.threshold FDR threshold for significance, if NULL return all hits
#' @param ... further arguments for \code{gprofiler2::gost()}
#' 
#' @author Alexander Toenges
#' 
#' @examples 
#' 
#' gosty(InputGenes = c("Ets1", "Cebpa", "Cd14", "Gapdh", "Spi1", "Cebpe"))
#' 
#' @export
gosty <- function(InputGenes, 
                  BackgroundGenes = NULL, 
                  Species = "mmusculus",
                  Sources = c("KEGG", "REAC"),
                  FDR.threshold = 0.05,
                  ...){
  
  gg <- gost(query = InputGenes, 
             custom_bg = BackgroundGenes,
             organism = Species, 
             exclude_iea = TRUE,
             evcodes = TRUE,
             user_threshold = FDR.threshold,
             sources = Sources, 
             ...)$result
  
  if(is.null(gg)) return(NULL)
  
  gg <- data.frame(Term = gg$term_name, 
                   pvalue = gg$p_value, 
                   isize = gg$intersection_size, 
                   tsize = gg$term_size, 
                   Source = gg$source, 
                   Genes = gg$intersection) 
  
  #/ Sort genes alphabetically:
  gg$Genes <- unlist(lapply(gg$Genes, function(x) paste(sort(strsplit(x, split = ",")[[1]]), collapse = ",")))
  
  #/ Sort by pvalue:
  gg <- gg%>%arrange(pvalue)
  
  return(gg)
  
}

