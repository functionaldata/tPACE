#' Normalise a curve to a particular area, by multiplication with a factor
#' 
#' Normalise a curve such that \\int{yNew}dx = area (according to trapezoid integration)
#' 
#' @param y values of curve at time-points x
#' @param x design time-points (default: seq(0,1, length.out=length(y)))
#' @param area value to normalise the curve onto (default: 1)
#'
#' @return yNew values of curve at times x such that [\\int{yNew}dx = area]
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

