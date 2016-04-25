# #' Check if NaN are present in the data and if yes remove them
# #' 
# #' Check if there are problems cause by missing values with the form and basic structure of the functional data 'Ly' and the recorded times 'Lt'.
# #' 
# #' @param Ly is a n-by-1 list of vectors
# #' @param Lt is a n-by-1 list of vectors

HandleNumericsAndNAN <- function(Ly,Lt){
 
  # Check for the presense of NA and remove them (if they exist) from the two lists in a pairwise manner
  if( any(is.na(unlist(Lt))) ||  any(is.na(unlist(Ly))) ){
   
    helperF <- function(x) which(!is.na(unlist(x)))
    L <- list(); for(j in 1:length(Ly)) L[[j]] = list(Ly[[j]],Lt[[j]])
    validIndexes = lapply(L, function(x) intersect(helperF(x[1]), helperF(x[2]) ))

    Ly = lapply(1:length(Ly), function(i) Ly[[i]][validIndexes[[i]]])
    Lt = lapply(1:length(Ly), function(i) Lt[[i]][validIndexes[[i]]])

    if( any(unlist(lapply(Ly, function(x) length(x) == 0))) ){
       stop('Subjects with only NA values are not allowed.\n')
    }
    
    ni_y = unlist(lapply(Ly,function(x) sum(!is.na(x))))
    if(all(ni_y == 1)){  
      stop("FPCA is aborted because the data do not contain repeated measurements after removing NA values.\n"); 
    }
  }
  

  
  # Force the data to be list of numeric members
  Ly <- lapply(Ly, as.numeric) 
  Lt <- lapply(Lt, as.numeric)
  Lt <- lapply(Lt, signif, 14)
  return( inputData <- list(Ly=Ly, Lt=Lt));

}
