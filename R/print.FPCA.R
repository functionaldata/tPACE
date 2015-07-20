print.FPCA <- function(obj){

  cat("Functional Principal Components Object for", tolower(obj$optns$dataType), "data.\n\n")
  cat("The optimal number of components selected is:", length(obj$lambda),"and \nthe first k (<=3) corresponding eigenvalues are: ");
  if ( length(obj$lambda) < 4) { 
    cat( round(obj$lambda,3) ,"\n");
  } else {
    cat( round(obj$lambda[1:3],3) ,"\n")
  }
}


 
