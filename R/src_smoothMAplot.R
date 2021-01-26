#' A function to produce MA-plots.
#' 
#' The function takes average log2-expression as X.values and logFC as Y.values and then uses
#' smoothScatter() to produce MA-plots. Due to smoothScatter() resulting pdfs are typically small
#' which is the reason we use this function rather than plot() or anything similar.
#' Optionally, if one provides values via P.values, data points whose P.values are below Sig.Thresh
#' will be colored in Sig.Color. 
#' 
#' @param X.value log2 average expression values
#' @param Y.value log2 fold changes
#' @param P.values optional p-values from differential analysis
#' @param Y.limits the limits of the y-axis, e.g. like c(-2,2), if NULL then automatically determined based on quantiles
#' to capture most data points while avoiding excessive axis limits due to outliers
#' @param Plot.xlab x-axis label
#' @param Plot.ylab y-axis label 
#' @param Plot.title title of the plot, 
#' @param Plot.title.size cexsize of title
#' @param Plot.cexlab cex size of the axis labels
#' @param Plot.cexaxis cex size of the axis marks
#' @param Plot.ColorRamp a color ramp for the plot, by default use the standard from smoothScatter() so white to blue
#' @param Points.cex cex size of the data points below the Sig.Thresh
#' @param Sig.Thresh threshold to select significant data points based on P.values
#' @param Sig.Color color to highlight significant data points, default is firebrick
#' @param Add.Legend logical, whether to add a legend at topright position with a summary of the number of data points 
#' above and below Sig.Thresh and the smallest arithmetic (not log2) fold change that was still significant at the given p-value cutoff
#' @param use.loessFit logical, whether to use loessFit for a fit between Yval and Xval
#' @param use.loessFit.Color color of that abline
#' @param abline0 logical, whether to add a y=0 abline
#' @param abline2 a numeric value, if set then will plot a horizontal abline at +/- that value
#' 
#' @details 
#' 
#' @author Alexander Toenges
#' 
#' @export
smoothMAplot  <-   function(X.value,
                            Y.value,
                            Y.limits = NULL,
                            y.quantile.low = 0.001,
                            y.quantile.top = 0.999,
                            Plot.xlab = "baseMean",
                            Plot.ylab = "logFC",
                            Plot.title = NULL,
                            Plot.title.size = 1,
                            Plot.cexlab = 1.5,
                            Plot.cexaxis = 1.5,
                            Plot.ColorRamp = NULL,
                            Points.cex = 0.2,
                            use.loessFit = FALSE,
                            use.loessFit.Color = "red",
                            abline0 = TRUE,
                            abline2 = NULL){
  
  #### Plot parameters: 
  x.value <- X.value
  y.value <- Y.value
  
  if(use.loessFit){
    
    # see section 4.6.3 of the csaw vignette, code inspired (taken) from there, thanks Aaron :)
    # https://www.bioconductor.org/packages/release/workflows/vignettes/csawUsersGuide/inst/doc/csaw.pdf
    lfit <- limma::loessFit(y = y.value, x = x.value, method = "lowess")
    o <- order(x.value)
    loess.toPlot <- data.frame(X=x.value[o], Y=lfit$fitted[o])
    
  }
  
  ## Definition of y-axis range, by default (if NULL) use quantiles 
  ## to avoid extreme axis ranges due to outliers:
  if (is.null(Y.limits)){
    Y.limits <- as.numeric( c( floor(quantile(y.value, y.quantile.low, na.rm=TRUE)), 
                               ceiling(quantile(y.value, y.quantile.top, na.rm=TRUE))))
  }
  
  ## Plot title:
  if(is.null(Plot.title)) title.basic <- "MAplot"
  if(!is.null(Plot.title)) title.basic <- Plot.title
  
  if(!is.null(Plot.ColorRamp)){
    cramp <- Plot.ColorRamp
  } else cramp <- colorRampPalette(c("white", blues9))
  
  ## the basic MA-plot:
  par(bty="L")
  smoothScatter(x        = x.value, 
                y        = y.value,
                xlab     = Plot.xlab,
                ylab     = Plot.ylab, 
                main     = title.basic,
                ylim     = Y.limits, 
                cex.main = Plot.title.size,
                cex.lab  = Plot.cexlab, 
                cex.axis = Plot.cexaxis,
                colramp  = cramp)
  
  ## add data points beyond Y.limitsits as trianges at the Y.limitsits
  lower.limit <- par("usr")[3]
  upper.limit <- par("usr")[4]
  
  my.low <- x.value[which(y.value < lower.limit)]
  my.high <- x.value[which(y.value > upper.limit)]
  
  if(length(my.low) > 0) {
    points(x = my.low, 
           y = rep( (lower.limit - lower.limit*0.025), length(which(y.value < lower.limit))),
           pch=17, cex = 0.8, col = cramp(1000)[1000]) 
  }
  
  if(length(my.high) > 0) {  
   points(x = my.high, 
           y = rep( (upper.limit - upper.limit*0.025), length(which(y.value > upper.limit))),
           pch=17, cex = 0.8, col=cramp(1000)[1000]) 
  }
  
  if(use.loessFit) lines(loess.toPlot$X,loess.toPlot$Y, col = use.loessFit.Color)
    
  if(abline0) abline(h=0, lty=2)
  if(!is.null(abline2)){
    abline(h=abline2, lty=2)
    abline(h=-abline2, lty=2)
  }
}
