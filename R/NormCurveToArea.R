#' Normalize a curve to a particular area, by multiplication with a factor
#' 
#' Normalize a curve such that its integration over the design time-points equals a particular value (according to trapezoid integration).
#' 
#' @param y values of curve at time-points \code{x}
#' @param x design time-points (default: \code{seq(0,1, length.out=length(y))})
#' @param area value to normalize the curve onto (default: 1)
#'
#' @return values of curve at times \code{x} such that its integration over \code{x} equals \code{area}.
#' @export

NormCurvToArea <- function(y, x = seq(0, 1, length.out = length(y)), area = 1){

  if( length(x) != length(y)){
    stop("'x' and 'y' must have the same length.")
  }                     
  if( length(y) < 2 ){
          stop("No area is defined for a single measurement.")
  }
  yNew = area * y / trapzRcpp(X = x, Y = y);
  return(yNew)
}

