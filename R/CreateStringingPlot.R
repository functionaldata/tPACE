#' Create plots for observed and stringed high dimensional data
#' 
#' The function produces the following three plots:
#' 1) A plot of predictors (standardized if specified so during stringing) in original order for a subset of observations;
#' 2) A plot of predictors in stringed order for the same subset of observations;
#' 3) A plot of the stringing function, which is the stringed order vs. the original order.
#' 
#' @param stringingObj A stringing object of class "Stringing", returned by the function Stringing.
#' @param subset A vector of indices or a logical vector for subsetting the observations. If missing,  first min(n,50) observations will be plotted where n is the sample size.
#' @param ... Other arguments passed into matplot for plotting options
#' @examples 
#' set.seed(1)
#' n <- 50
#' wiener = Wiener(n = n)[,-1]
#' p = ncol(wiener)
#' rdmorder = sample(size = p, x=1:p, replace = FALSE)
#' stringingfit = Stringing(X = wiener[,rdmorder], disOptns = "correlation")
#' diff_norev = sum(abs(rdmorder[stringingfit$StringingOrder] - 1:p))
#' diff_rev = sum(abs(rdmorder[stringingfit$StringingOrder] - p:1))
#' if(diff_rev <= diff_norev){
#'   stringingfit$StringingOrder = rev(stringingfit$StringingOrder)
#'   stringingfit$Ly = lapply(stringingfit$Ly, rev)
#' }
#' CreateStringingPlot(stringingfit, 1:20)
#' 
#' @export

CreateStringingPlot <- function(stringingObj, subset, ...){
  if(class(stringingObj) != "Stringing"){
    stop("Invalid input class for the stringing object. Need to be returned from the Stringing function.")
  }
  n <- nrow(stringingObj$Xin)
  p <- ncol(stringingObj$Xin)
  if(missing(subset)){
    subset <- seq_len(min(n, 50))
  }
  subset = subset[which(subset > 0)]
  subset = unique(subset)
  if( !all(subset %in% seq_len(n)) ){
    stop("Invalid subset of subjects specified for generating plots. Need to be a subset of the set of row indices.")
  }
  ylabel = "Standardized Predictors"
  X = stringingObj$Xstd
  if(is.null(stringingObj$Xstd)){
    X = stringingObj$Xin
    ylabel = "Observed Predictors"
  }
  
  inargs = list(...) # additional options for plotting
  
  defaultColPalette = palette()
  args1 <- list( xlab= 'Observed Order', ylab= ylabel, col = defaultColPalette, pch=1)    
  args1[names(inargs)] <- inargs
  
  args2 <- list( xlab= 'Stringed Order', ylab= ylabel, col = defaultColPalette, pch=1)    
  args2[names(inargs)] <- inargs
  
  plotx = matrix(rep(seq_len(p), length(subset)), nrow = length(subset), byrow = TRUE)
  ploty1 = X[subset, ,drop = FALSE]
  ploty2 = X[subset, stringingObj$StringingOrder, drop = FALSE]
  
  # 1. plot for the original observations
  # Make canvas
  do.call(plot, c(list(x=plotx, y=ploty1, type='n'), args1))
  # plot the observed/standardized predictors
  do.call(matplot, c(list(x=t(plotx), y=t(ploty1), type='l', add = TRUE), args1))
  readline("Press enter to continue.")
  
  # 2. plot for the stringed observations
  # Make canvas
  do.call(plot, c(list(x=plotx, y=ploty2, type='n'), args2))
  # plot the stringed predictors
  do.call(matplot, c(list(x=t(plotx), y=t(ploty2), type='l', add = TRUE), args2))
  readline("Press enter to continue.")

  # 3. plot the stringing function
  plot(seq_len(ncol(X)), stringingObj$StringingOrder, 
       xlab="Observed Order", ylab="Stringed Order", pch = 18)
  par(mfrow=c(1,1))

  invisible()
}
