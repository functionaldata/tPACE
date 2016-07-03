#' Print an FSVD object
#'
#' Print a simple description of an FSVD object
#'
#' @param x An FSVD object.
#' @param ... Not used.
#'
#' @method print FSVD
#' @export
print.FSVD <- function(x, ...){
  obj = x;
  thisDataType <- NULL
  if(obj$optns$SVDopts$dataType1 == 'Dense' && obj$optns$SVDopts$dataType2 == 'Dense'){
    thisDataType <- 'Dense'
  } else {
    thisDataType <- 'Sparse'
  }
  if(obj$optns$SVDopts$dataType1 == 'DenseWithMV' && obj$optns$SVDopts$dataType2 == 'DenseWithMV'){
    thisDataType <- 'DenseWithMV'
  }
    
  message(paste0("Functional Singular Value Decomposition object for", tolower(thisDataType), "data.\n\n"))
  message(paste0("The optimal number of components selected is:", length(obj$sValues),"and \nthe first K (<=3) singular values are: "))
  if ( length(obj$sValues) < 4) { 
    message(paste0( round(obj$sValues,3) ,"\n"))
  } else {
    message(paste0(round(obj$sValues[1:3],3) ,"\n"))
  }
}


 
