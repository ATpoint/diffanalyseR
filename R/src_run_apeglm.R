#' A function to test contrasts via DESeq2/lfcShrink/apeglm
#' 
#' The function is designed to use apeglm for shrinking LFCs while automatically 
#' releveling the dds object to make the relevant coef accessable via resultsNames
#' so apeglm can access the MLEs, see details.
#' 
#' @param dds a DESeqDataSet() that contains size factors and dispersion estimates, see examples
#' @param contrasts named lists of contrasts in DESeq2 style, see details
#' @param lfc lfc threshold to use as Null
#' @param svalue logical, whether to return s-values rather than pvalues/padj. Will be TRUE automatically if lfc argument is not zero.
#' The latter is hardcoded in \code{DESeq2::lfcShrink}. 
#' @param workers number of BiocParallel workers
#' @param verbose logical, whether to allow status messages
#' 
#' @author Alexander Toenges
#' 
#' @details 
#' By design the \code{DESeq2::lfcShrink(type="apeglm")} function uses the coef argument rather than contrasts. 
#' Therefore the \code{resultsNames(dds)} must contain the respective coef so apeglm can access the MLEs of the effect sizes.  
#' This function here will relevel the factor levels of the dds object to make all contrasts accessable via coef.
#' A more illustrative example is provided in <https://www.biostars.org/p/448959/#484944>.
#' 
#' Say one makes a design as \code{~condition} with three levels A/B/C then the expected contrasts
#' argument of this function would be a named list like...
#' 
#' $A_vs_B
#' [1] "condition" "A" "B"
#' $A_vs_C
#' [1] "condition" "A" "C"
#' $B_vs_C
#' [1] "condition" "B" "C"
#' 
#' ...which is the DESeq2 style of writing down contrasts.
#' This function will then loop through the contrasts, releveling the factor levels automatically (if necessary)
#' to allow apeglm shrinkage for all comparisons. After a releveling the Wald test must be re-run, but not of the dispersion estimation because
#' the design is the same. 
#' 
#' So far this has only been tested for simple designs such as the one above. 
#' Not tested for more complex designs such as interactions etc.
#' If one uses the \code{lfc} argument to test against a Null other than zero then s-values will always be returned rather than p/padj.
#' This is hardcoded in DESeq2 and cannot be changed. The vignette of DESeq2 recommends a smaller cutoff for s- rather the padj, e.g. 0.005.
#' Final output is a named list with one data.frame per contrast storing the DE statistics.
#' 
#' @examples
#' # some dummy data:
#' dds <- DESeq2::makeExampleDESeqDataSet()
#' dds$condition <- factor(unlist(lapply(seq(1,3),function(x) rep(LETTERS[x], 4))))
#' contrasts <- diffanalyseR::makeAllContrasts(levels(dds$condition), deseq2 = TRUE, name = "condition")
#' dds <- DESeq2::estimateSizeFactors(dds) %>% DESeq2::estimateDispersions(.)
#' res <- run_apeglm(dds,contrasts, lfc=.5)
#' 
#' @export
run_apeglm <- function(dds, contrasts, lfc=.5, 
                       svalue=TRUE, workers=1, verbose=TRUE)
  
{
  
  ########################################
  # Checks
  ########################################
  
  pkg <- c("DESeq2", "apeglm") 
  if(sum(pkg %in% rownames(installed.packages())) != length(pkg)) 
    stop("Make sure both DESeq2 and apeglm are installed!", call. = FALSE)
  
  if(is.null(DESeq2::sizeFactors(dds)) & ! "normalizationFactors" %in% SummarizedExperiment::assayNames(dds))
    stop("dds contains neither size- nor normalization factors", call. = FALSE)
  if(is.null(mcols(dds)$dispersion)) stop("dds contains no dispersion estimates", call. = FALSE)
  invisible(match.arg(arg = class(contrasts)[1], choices = c("list")))
  invisible(match.arg(arg = class(lfc), choices = c("numeric", "integer")))
  
  ########################################
  # Main
  ########################################
  
  #/ data.frame summarizing the contrasts:
  df <- data.frame(t(data.frame(contrasts)))
  df$X4 <- paste0(df[,1], "_", df[,2], "_vs_", df[,3])
  df <- df %>% arrange(X3)
  nm <- unique(df[,1])
  
  #/ BPPARAM:
  if(workers>1) {
    do.parallel=TRUE 
    BP=BiocParallel::MulticoreParam(workers)
  } else {
    do.parallel=FALSE
    BP=BiocParallel::SerialParam()
  }
  
  #/ loop through the contrasts:
  res <- list()
  for(i in 1:nrow(df)){
    
    dds[[nm]] <- relevel(dds[[nm]], df[i,3])
    
    #/ Run Wald if necessary for that comparison:
    if(length(DESeq2::resultsNames(dds))==0 | !df[i,4] %in% DESeq2::resultsNames(dds)){
      if(verbose) message("Releveling to ", df[i,3])
      dds <- DESeq2::nbinomWaldTest(object = dds, quiet = TRUE)
    }
    
    # lfcShrink:
    if(verbose) message("=> Testing ", df[i,4])
    tt <- suppressMessages(DESeq2::lfcShrink(dds = dds, coef = which(DESeq2::resultsNames(dds)==df[i,4]), 
                            lfcThreshold = lfc, quiet = !verbose, svalue = svalue, type = "apeglm",
                            parallel = do.parallel, BPPARAM = BP)) %>% 
            data.frame(.) %>% 
            data.frame(Gene=rownames(.), .) %>% 
            na.omit
    
    tt$baseMean <- log2(tt$baseMean+1)
    
    if(svalue){
      colnames(tt) <- c("Gene","baseMean","logFC","lfcSE","svalue")
      tt <- tt[,c(1,3,2,4,5)]
    } else{
      tt <- tt[,c(1,3,2,4,5,6)]        
    }
    
    res[[i]] <- tt
    rm(tt)
    
  };  names(res) <- rownames(df)
  
  return(res)
  
}
