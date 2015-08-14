#' Check data format
#' 
#' Check the form and basic structure of the functional data 'y' and the recorded times 'tt'.
#' 
#' @param y is a n-by-1 list of vectors
#' @param tt is a n-by-1 list of vectors
#' @return logical
#' @examples 
#' 1 + 3
#' @export

CheckData = function(y,t){
  
  if(!is.list(y)){
    cat('Error: y should be list \n')
    return(TRUE);
  }
  if(!is.list(t)){
    cat('Error: t should be list \n')
    return(TRUE);
  }
  if(any(is.na(unlist(y)))){
    cat('Error: y cannot contain NA/NaN entries\n')
    return(TRUE);
  }
  if(any(is.na(unlist(t)))){
    cat('Error: t cannot contain NA/NaN entries\n')
    return(TRUE);
  } 
  ni_y = unlist(lapply(y,function(x) sum(!is.na(x))))
  if(all(ni_y == 1)){  
    cat("Error:FPCA is aborted because the data do not contain repeated measurements in y!\n"); 
    return(TRUE);    
  }
  ni_tt = unlist(lapply(t,function(x) sum(!is.na(x))))
  if(all(ni_tt == 1)){  
    cat("Error:FPCA is aborted because the data do not contain repeated measurements in t!\n"); 
    return(TRUE);    
  }   
  if( !all(unlist(lapply(y,function(x) class(x) %in% c('integer', 'numeric') ) ) ) ){
        cat("Error:FPCA is aborted because 'y' members are not all of class numeric! Try  \"lapply(y,function(x) class(x))\" to see the current classes. \n");     return(TRUE);
  }
 if( !all(unlist(lapply(t,function(x) class(x) %in% c('integer', 'numeric'))) ) ){
        cat("Error:FPCA is aborted because 't' members are not all of class numeric! Try  \"lapply(t,function(x) class(x))\" to see the current classes. \n");     return(TRUE);
  }

 if(any( unlist( lapply(t, function(x) length(x) != length(unique(x))))) ){
        cat("Error:FPCA is aborted because within-subject 't' members have duplicated values.  Try  \"which( unlist( lapply(t, function(x) length(x) != length(unique(x)))))\" to see potentially problematic entries. \n");     return(TRUE);
  }

  return(FALSE);
}

