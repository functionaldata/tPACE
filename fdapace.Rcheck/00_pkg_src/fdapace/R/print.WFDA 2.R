#' Print a WFDA object                                                                                                           
#'
#' Print a simple description of a WFDA object
#'
#' @param x A WFDA object.
#' @param ... Not used.
#'
#' @method print WFDA
#' @export
print.WFDA <- function(x, ...){
    obj = x;
    
    cat(paste0("Warped Functional Data Analysis object for ", length(obj$costs), " curves.\n\n"))
    cat(paste0("The penalty parameter used was: ", signif(obj$lambda,6), ", the warping functions ", 
               ifelse(obj$optns$isPWL,"are ", "are not"), "piece-wise linear \n",
               "and the pairwise warping was done using the ", obj$optns$choice, " averages of the warped curves.\n" ))
}

