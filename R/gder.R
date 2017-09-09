#' Derivative of the link function
#' 
#' Get the derivative of the link function
#' 
#' @param muy is a n-by-1 vector
#' @param family is a string specifying the distribution family of the variable; has to be 'poisson', 'binomial' or 'gaussian' (default: 'poisson')
#' @export

gder <- function(muy, family = 'poisson'){
  if(family=='poisson'){
    gder = exp(muy)
  } else{
    if( family=='binomial'){ 
      gder = exp(muy)/(1+exp(muy))^2
    } else {
      if(family == 'gaussian'){
        gder = muy
      } else {
        stop('family argument is wrong.')
      }
    }
  }
  return(gder)
}