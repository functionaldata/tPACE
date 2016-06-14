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
    
  cat("Functional Singular Value Decomposition object for", tolower(thisDataType), "data.\n\n")
  cat("The optimal number of components selected is:", length(obj$sValues),"and \nthe first k (<=3) singular values are: ");
  if ( length(obj$sValues) < 4) { 
    cat( round(obj$sValues,3) ,"\n");
  } else {
    cat( round(obj$sValues[1:3],3) ,"\n")
  }
}


 
