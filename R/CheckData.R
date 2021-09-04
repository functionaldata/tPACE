#' Check data format
#' 
#' Check if there are problems with the form and basic structure of the functional data 'y' and the recorded times 't'.
#' 
#' @param y is a n-by-1 list of vectors
#' @param t is a n-by-1 list of vectors
#' @export


CheckData = function(y,t){
  
  if(!is.list(y)){
    stop('y should be list \n')
  }
  if(!is.list(t)){
    stop('t should be list \n')
  }
  
  if( length(t) != length(y)){
    stop('t and y should have the same length \n')
  }
  
  ni_y = unlist(lapply(y,function(x) sum(!is.na(x))))
  if(all(ni_y == 1)){  
    stop("FPCA is aborted because the data do not contain repeated measurements in y!\n"); 
  }
  ni_tt = unlist(lapply(t,function(x) sum(!is.na(x))))
  if(all(ni_tt == 1)){  
    stop("FPCA is aborted because the data do not contain repeated measurements in t!\n"); 
  }   
  if( !all(unlist(lapply(y,function(x) typeof(x) %in% c('integer', 'double') ) ) ) ){
    stop("FPCA is aborted because 'y' members are not all of type double or integer! Try  \"lapply(y,function(x) typeof(x))\" to see the current types \n");
  }
  if( !all(unlist(lapply(t,function(x) typeof(x) %in% c('integer', 'double'))) ) ){
    stop("FPCA is aborted because 't' members are not all of type double or integer! Try  \"lapply(t,function(x) typeof(x))\" to see the current types \n");
  }
  
  if(any( unlist( lapply(t, function(x) length(x) != length(unique(x))))) ){
    stop("FPCA is aborted because within-subject 't' members have duplicated values.  Try  \"which( unlist( lapply(t, function(x) length(x) != length(unique(x)))))\" to see potentially problematic entries. \n");
  }
  if( any(sapply(t[seq_len(min(1001, length(t)))], is.unsorted, na.rm=TRUE)) ) {
    stop('Each vector in t should be in ascending order')
  }
  if(min(unlist(y),na.rm=TRUE)==-Inf){
    stop('There are entries in Ly which are -Inf')
  }
  if(max(unlist(y),na.rm=TRUE)==Inf){
    stop('There are entries in Ly which are Inf')
  }
}

