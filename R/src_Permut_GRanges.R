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
                           return.values=TRUE){
  
  invisible(match.arg(class(featureset), "GRanges"))
  invisible(match.arg(class(foreground), "GRanges"))
  invisible(match.arg(class(background), "GRanges"))
  invisible(match.arg(class(restrictset), c("GRanges", "NULL")))
  
  #/ restrictset is provided. give each entry a unique name and use this
  #/ as seqnames
  if(!is.null(restrictset)){
    
    restrictset$ID <- paste0("ID", 1:length(restrictset))
    
    #/ take the center of the feature/fore/background and intersect with restrictset:
    makenewset <- function(oset, restr){
      olap <- suppressWarnings(findOverlaps(resize(oset,fix="center",width=1), restr))
      suppressWarnings(GRanges(seqnames = restr[olap@to]$ID, ranges = ranges(oset[olap@from])))
    }
    
    featureset=makenewset(oset=featureset, restr=restrictset)
    foreground=makenewset(oset=foreground, restr=restrictset)
    background=makenewset(oset=background, restr=restrictset)
    
  }
  
  # Permutate the background based on foreground length, get idx:
  idx <- lapply(1:nperm, function(x){
    sample(1:length(background), length(foreground), replace=FALSE)
  })
  
  #/ Sample from background based on the idx
  grl <- background[unlist(idx)] %>%
    plyranges::mutate(permut=rep(paste0("perm",1:length(idx)),each=length(foreground)))
  
  #/ find overlaps between featureset and permutated background (=grl)
  s.random <- suppressWarnings(findOverlaps(featureset, grl)) %>%
    data.frame(featureset[.@from], grl[.@to]) %>%
    select(c("seqnames","start","end","permut")) %>%
    distinct %>%
    pull(permut) %>%
    table %>%
    as.numeric
  
  #/ in case some overlaps were empty add zeros:
  s.random <- c(s.random, rep(0, nperm-length(s.random)))
  
  #/ Observed overlap:
  olap.fg <- findOverlaps(featureset, foreground)
  s.is    <- length(unique(olap.fg@from))
  per_element.foreground <- 
    table(olap.fg@to) %>%
    as.data.frame %>%
    merge(x=data.frame(A=seq(1,length(foreground))), y=.,
          by.x="A",by.y="Var1", all.x=TRUE) %>%
    mutate(Freq=ifelse(is.na(Freq), 0, Freq)) %>%
    pull(Freq)
  
  olap.bg <- findOverlaps(featureset, background)
  s.bg    <- length(unique(olap.bg@from))
  per_element.background <- 
    table(olap.bg@to) %>%
    as.data.frame %>%
    merge(x=data.frame(A=seq(1,length(foreground))), y=.,
          by.x="A",by.y="Var1", all.x=TRUE) %>%
    mutate(Freq=ifelse(is.na(Freq), 0, Freq)) %>%
    pull(Freq)
  
  #/ Summary per foreground region:
  
  #/ construct output list:
  l <- list(pval_greater=(1+sum(s.random >= s.is))/(nperm+1), 
            pval_less=(1+sum(s.random <= s.is))/(nperm+1), 
            observed=s.is, 
            CI=tryCatch(expr=as.numeric(t.test(s.random,conf.level = conf.level)$conf.in),
                        error=function(x) return(NA)),
            individual_foreground=per_element.foreground)
  
  if(return.values) l[["s.random"]] <- s.random
  
  return(l)
  
}
