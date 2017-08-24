#' Inverse of the link function
#' 
#' Get the inverse of the link function
#' 
#' @param muy is a n-by-1 vector
#' @param family is a string specifying the distribution family of the variable; has to be 'poisson', 'binomial' or 'gaussian' (default: 'poisson')
#' @export

ginver <- function(muy, family = 'poisson'){
  if(family=='poisson'){
    muy[muy<=0] = exp(-5)
    ginver = log(muy)
  } else{
    if( family=='binomial'){
      muy[muy<=0] = exp(-5);
      muy[muy>=1] = 1-exp(-5);
      ginver = log(muy/(1-muy));
    } else {
      if(family == 'gaussian'){
        ginver = muy
      } else {
        stop('family argument is wrong.')
      }
    }
  }
  return(ginver)
}
