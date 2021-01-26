#' Convert strsplit to data.frame
#' 
#' Assuming coordinates such as chr:start-end, function splits at delimiters and produces
#' data.frame output. Uses `data.table::tstrsplit` internally.
#' 
#' @param x input string(s)
#' @param pattern pattern to split string at. Can be regexed.
#' @param chr.name name for chr
#' @param start.name name for start
#' @param end.name name for end
#' 
#' @details 
#' 
#' @author Alexander Toenges
#' 
#' @examples 
#' strs <- c("chr1:1-10", "chr2:2-20")
#' string2chr(strs, ":|-")
#' 
#' @export
string2chr <- function(x,
                         pattern = ":|-",
                         chr.name = "chr",
                         start.name = "start",
                         end.name = "end"){
  
  splitted <- data.table::tstrsplit(x, pattern)
  if(length(splitted)!=3) stop("Split did not produce three elements")
  df<-data.frame(A=splitted[[1]], B=splitted[[2]], C=splitted[[3]])
  colnames(df) <- c(chr.name, start.name, end.name)
  df
  
}
 