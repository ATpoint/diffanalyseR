#' Combine AveTxLen with norm.factors into offsets
#' 
#' Wrapper around https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#edger
#' Combine edgeR normalization factors with average transcript length and obtain norm. offsets for downstream use.
#' 
#' @param cts count matrix
#' @param len length matrix
#' @param norm_method norm method from edgeR::calcNormFactors
#' 
#' @details 
#' 
#' @author Alexander Toenges
#' @export
NormOffsetsFromTxLength <- function(cts, len, norm_method = "TMM"){
  
  normMat <- len
  
  # Obtaining per-observation scaling factors for length, adjusted to avoid
  # changing the magnitude of the counts.
  normMat <- normMat/exp(rowMeans(log(normMat)))
  normCts <- cts/normMat
  
  # Computing effective library sizes from scaled counts, to account for
  # composition biases between samples.
  eff.lib <- edgeR::calcNormFactors(normCts, method = norm_method) * colSums(normCts)
  
  # Combining effective library sizes with the length factors, and calculating
  # offsets for a log-link GLM.
  normMat <- base::sweep(normMat, 2, eff.lib, "*")
  normMat <- log(normMat)
  
  # Calculate offsets and add to the SE.
  return(edgeR::scaleOffset(DGEList(counts = cts), normMat)$offset)
  
}