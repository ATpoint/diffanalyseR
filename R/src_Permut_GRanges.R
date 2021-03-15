#' A function to generate permutation-based p-values based on overlap between GRanges objects
#'
#' Output is a named list with p-values, observed overlap and confidence intervals for 
#' the permutation-based overlap test.
#'
#' @param featureset GRanges object to overlap with the below objects, see details.
#' @param foreground GRanges object containing regions that will be used to generate 
#' the observed overlaps with the featureset set
#' @param background GRanges object containing regions that will be used to generate 
#' the permutated/random overlaps with the featureset set
#' @param restrictset GRanges to restrict overlapping operations between featureset and fore/background to.
#' This could e.g. be TAD boundaries and only overlaps between featureset and fore/background within the same TAD would be counted.
#' Features outside TADs are ignored.
#' @param nperm number of permutations. Smallest possible p-value is \code{1/(nperm+1)}
#' @param return.values logical, whether to return the number of overlaps between permuted backgrounds and featureset
#' as part of the output list. Can be useful for plotting/visualization, is not ordered though.
#' @param conf.level the confidence level to calculate 
#' 
#' @details 
#' Underlying question is whether the overlap for a given set of features (featureset), e.g.
#' ChIP-seq peaks is greater between a given foreground set, e.g. certain TSS, compared to a 
#' background, e.g. all annotated TSS. The permutation-based test performs overlap counting between
#' the foreground and permutated background (number of permutations by default 5000). Pvalues are then calculated
#' which test whether the background overlaps are equal or more extreme than the foreground overlaps.
#' The smallest possible pvalue is defined by 1/(nperm+1).
#' Optionally one can provide a restrictset, e.g. TAD coordinates, so only intersections between features
#' are counted that overlap the same restrictset interval.
#' 
#' @examples 
#' library(airway)
#' data(airway)
#' 
#' background=reduce(unlist(rowRanges(airway)[c(1:500)]))
#' set.seed(2021)
#' foreground=background[sample(1:length(background), 50, replace=FALSE)]
#' featureset=foreground
#' Permut_GRanges(featureset=featureset,foreground=foreground,background=background)
#' 
#' @author Alexander Toenges
#' 
#' @export
Permut_GRanges <- function(featureset, foreground,  
                           background, restrictset=NULL,
                           nperm=5000, conf.level=.95,
                           return.values=FALSE)
{
  
  invisible(match.arg(class(featureset), "GRanges"))
  invisible(match.arg(class(foreground), "GRanges"))
  invisible(match.arg(class(background), "GRanges"))
  invisible(match.arg(class(restrictset), c("GRanges", "NULL")))
  
  #/ If a restrictset is provided make a little trick to speed up computations,
  #/ simply give each element of restrictset an "ID" and use this as the
  #/ chromosome identifier. That way we do not need to run the intersection
  #/ for each element separately.
  if(!is.null(restrictset)){
    
    restrictset$ID <- paste0("ID", 1:length(restrictset))
    makenewset <- function(oset, restr){
      olap <- findOverlaps(oset, restr)
      suppressWarnings(GRanges(seqnames = restr[olap@to]$ID, ranges = ranges(oset[olap@from])))
    }
    featureset=makenewset(oset=featureset, restr=restrictset)
    foreground=makenewset(oset=foreground, restr=restrictset)
    background=makenewset(oset=background, restr=restrictset)
    
  }
  
  # Get idx for the permutations without replacement
  idx <- lapply(1:nperm, function(x){
    sample(1:length(background), length(foreground), replace=FALSE)
  })
  
  #/ Sample from background based on the idx
  grl <- background[unlist(idx)]
  
  #/ add a unique name to each set of permutations
  grl$permut <- unlist(lapply(1:nperm, function(x) rep(paste0("perm",x), length(foreground))))
  
  #/ find overlaps between foreground and permutated background (=grl)
  fo <- suppressWarnings(findOverlaps(featureset, grl))
  
  s.random <- as.numeric(table(data.frame(featureset[fo@from], grl[fo@to])$permut))
  
  #/ if s.random is smaller than nperm means there were zero overlaps for some
  #/ permutations, fill them up with zeros
  s.random <- c(s.random, rep(0, nperm-length(s.random)))
  
  #/ the actual (=observed) overlap:
  s.is <- length(findOverlaps(featureset, foreground)@from)
  
  #/ get permutation-based p-values testing the Null that the permutations return
  #/ the same or more extreme (or less extreme) values than the observed overlap:
  pval_greater <- (1+sum(s.random >= s.is))/(nperm+1)
  pval_less    <- (1+sum(s.random <= s.is))/(nperm+1)
  
  #/ CI:
  ci <- tryCatch(expr=as.numeric(t.test(s.random,conf.level = conf.level)$conf.in),
                 error=function(x) return(NA))
  ci[ci<0] <- 0
  
  l <- list(pval_greater=pval_greater, pval_less=pval_less, observed=s.is, CI=ci)
  if(return.values) l[["s.random"]] <- s.random
  return(l)
  
}