#' Row scaling with matrixStats
#' 
#' An efficient function for row scaling using matrixStats, see reference for origin.
#'
#' @param x input matrix
#' @param center logical, whether to center
#' @param scale logical, whether to scale
#' @param add_attr logical, whether to add scale/center the attributes
#' @param rows rows to subset for
#' @param cols cols to subset for
#' 
#' @details
#' 
#' @author Alexander Toenges
#' 
#' @examples 
#' 
#' @references 
#' A Faster Scale Function (2016) https://www.r-bloggers.com/a-faster-scale-function/
#'
#' @export
rowScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = FALSE,
                    rows = NULL,
                    cols = NULL) {
  
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  
  ################
  # Get the column means
  ################
  cm = rowMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = matrixStats::rowSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}