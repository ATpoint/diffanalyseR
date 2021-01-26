#' Automate the glmQLFit/glmTreat functions from edgeR
#' 
#' Function assumes a DGEGLM object from glmQLFit and will then produce DE results based on a coef or contrast.
#' 
#' @param fit DGEGLM, output of glmQLFit
#' @param contrast a numeric contrast so a column from edgeR::makeContrasts()
#' @param coef a coefficient to test, ignored it contrast is set
#' @param lfc a lfc threshold for glmTreat, if 0 falls back to glmQLTest
#' @param ID the name of the column that will store the genes/regions/identifiers for each row
#'
#' @author Alexander Toenges
#' 
#' @examples 
#' nlibs <- 4
#' ngenes <- 1000
#' dispersion.true <- 1/rchisq(ngenes, df=10)
#' design <- model.matrix(~factor(c(1,1,2,2)))
#' y <- rnbinom(ngenes*nlibs,mu=20,size=1/dispersion.true)
#' y <- matrix(y,ngenes,nlibs)
#' rownames(y) <- paste0("Gene", seq(1,nrow(y)))
#' d <- DGEList(y)
#' d <- calcNormFactors(d)
#' d <- estimateDisp(d, design)
#' fit <- glmQLFit(d, design)
#' head(res_edgeR(fit=fit,coef=2))
#' 
#' @export
res_edgeR <- function(fit, contrast = NULL, coef = NULL, lfc = log2(1.2), ID = "Gene"){
  
  invisible(match.arg(class(fit), "DGEGLM"))
  invisible(match.arg(class(contrast), c("numeric", "NULL")))
  invisible(match.arg(class(coef), c("numeric", "NULL")))
            
  if(!is.null(contrast) & sum(contrast == 0) == length(contrast)) stop("contrast is all-zero", call. = FALSE)
  
  # use glmTreat even if lfc=0, will fall back to glmQLTest
  sink(file="/dev/null")
  qlf <- edgeR::glmTreat(glmfit = fit, contrast = contrast, coef = coef, lfc = lfc)
  sink()
  
  tt  <- edgeR::topTags(qlf, n=Inf, adjust.method="BH", sort.by="none")$table %>% data.frame(ID=rownames(.), .)
  colnames(tt)[colnames(tt)=="ID"] <- ID
  colnames(tt)[colnames(tt)=="logCPM"] <- "baseMean"
  colnames(tt)[colnames(tt)=="PValue"] <- "pvalue"
  
  return(tt)
    
}
