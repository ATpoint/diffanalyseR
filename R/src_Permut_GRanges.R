#' A function to generate permutation-based p-values based on overlap between GRanges objects
#'
#' Output is a named list with p-values, observed overlap and confidence intervals.
#'
#' @param featureset GRanges object to overlap with the below objects, see details.
#' @param foreground GRanges object containing regions that will be used to generate 
#' the observed overlaps with the featureset set
#' @param background GRanges object containing regions that will be used to generate 
#' the permutated/random overlaps with the featureset set
#' @param restrictset GRanges that are used to restrict overlapping operations between featureset and fore/background to.
#' This could e.g. be TAD boundaries and only overlaps between featureset and fore/background within the same TAD would be counted.
#' Features outside TADs are currently ignored.
#' @param nperm number of permutations. Smallest p-value is calculated by 1/(nperm+1)
#' @param return.bg_overlaps logical, whether to return the number of overlaps between permuted backgrounds and featureset
#' as part of the output list. Can be useful for plotting/visualization.
#' @param conf.level the confidence level to calculate 
#' @param cores number of mclapply workers for the permutations
#' 
#' @details 
#' Test whether a set of genomic regions, here called foreground, is statistically associated with a set
#' of genomic features, here called featureset, e.g. foreground being a set of TSS, and featureset being ChIP-seq
#' TF peaks. The question could e.g. be whether these TSS and its vicinity (say 50kb) show more overlap with the TF peaks
#' than a background which could be all annotated TSS or all expressed TSS in this experiment. For this we would provide 
#' as foreground the TSS coordinates, resizes by 25kb in each direction. Based on these we count the observed overlap between 
#' featureset and foreground as the number of unique overlapping elements of featureset with foreground.
#' We compare this with the expected overlap based on the background. The background could e.g. be all annotated TSS 
#' (or all expressed TSS) extended by the same window size as foreground. Foreground can (and probably should) be part of background.
#' For the expected overlap based on background we randomly draw (without replacement) from the background the same number of entries as 
#' there are in foreground and then do the same overlap counting as above. We repeat this nperm times to define the expected
#' overlap and calculate confidence intervals. p-values are then calculated based on the formula:  
#' \code{(1+sum(s.random >= s.is))/(nperm+1)} where \code{s.random} are the overlaps between background and featureset,
#' \code{s.is} is the observed overlap between foreground and featureset and nperm is the number of permutations.
#' The returned `pval_greater` is then the probability for the Null of foreground not being enriched for overlaps.
#' The `pval_less` is the probability for the Null that the foreground is not depleted for overlaps. The formula above
#' is then `<=` rather than `>=`. The returned confidence intervals indicate the bounds of the expected overlaps based on the
#' background set. Output is returned as a named list.
#' 
#' @examples 
#' # A simple example, using airway to get some coordinates:
#' library(airway)
#' data(airway)
#' # simply use the first 500 genes to define a background set
#' background=reduce(unlist(rowRanges(airway)[c(1:500)]))
#' # as a primitive example we define the first 50 as the foreground
#' set.seed(2021)
#' foreground=background[sample(1:length(background), 50, replace=FALSE)]
#' # and as features we also use the first 50 so obviously these will be significantly
#' # associated with foreground as they're the same, this is just for illustration obviously
#' featureset=foreground
#' # unsurprisingly this is significant and the smallest possible pvalue (1/(nperm+1))
#' Permut_GRanges(featureset=featureset,foreground=foreground,background=background,nperm=100)
#' 
#' @author Alexander Toenges
#' 
#' @export
Permut_GRanges <- function(featureset,  
                           foreground,  
                           background,  
                           restrictset=NULL,
                           nperm=500,     
                           return.bg_overlaps=FALSE,
                           conf.level=.95,
                           cores = detectCores()){
  
  invisible(match.arg(class(featureset), "GRanges"))
  invisible(match.arg(class(foreground), "GRanges"))
  invisible(match.arg(class(background), "GRanges"))
  invisible(match.arg(class(restrictset), c("GRanges", "NULL")))
  
  if(!is.null(restrictset)){
    
    restrictset$ID <- paste0("ID", 1:length(restrictset))
    makenewset <- function(oset, restr){
      olap <- findOverlaps(oset, restr)
      GRanges(seqnames = restr[olap@to]$ID, ranges = ranges(oset[olap@from]))
    }
    featureset=makenewset(oset=featureset, restr=restrictset)
    foreground=makenewset(oset=foreground, restr=restrictset)
    background=makenewset(oset=background, restr=restrictset)
    
  }
  
  #/ Observed overlap so our Null
  s.is <- length(suppressWarnings(subsetByOverlaps(featureset, foreground)))
  
  #/ Background overlap function. Sample as many regions from background as there are foreground regions.
  .RandomFunction <- function(len.is, len.random, replace=FALSE) sample(1:len.random, len.is, replace = replace)
  
  len.is=length(foreground)
  len.random=length(background)
  
  if(len.is>len.random) warning("Foreground has more elements than background")
  
  #/ Count the overlaps between the featureset and the randomly-drawn elements from the background:
  s.random <- unlist(mclapply(seq(1,nperm), mc.cores = cores, 
                              function(x) length(suppressWarnings(subsetByOverlaps(featureset, 
                                                                                   background[.RandomFunction(len.is,len.random)])))))
  
  #/ get p-values:
  pval_greater <- (1+sum(s.random >= s.is))/(nperm+1)
  pval_less    <- (1+sum(s.random <= s.is))/(nperm+1)
  
  #/ CI:
  ci <- t.test(s.random,conf.level = conf.level)$conf.in
  
  l <- list(pval_greater=pval_greater, pval_less=pval_less, observed=s.is, CI=ci)
  if(return.bg_overlaps) l[["s.random"]] <- s.random
  return(l)
  
}