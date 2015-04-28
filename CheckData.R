CheckData = function(y,tt){
  
  # Check the form and basic structure of the functional data 'y'
  # and the recorded times 'tt'. 
  # y : n-by-1 list of vectors
  # tt : n-by-1 list of vectors 
  
  if(!is.list(y)){
    cat('Error: y should be list \n')
    return(TRUE);
  }
  if(!is.list(tt)){
    cat('Error: tt should be list \n')
    return(TRUE);
  }
  if(any(is.na(unlist(y)))){
    cat('Error: y cannont contain NA/NaN entries\n')
    return(TRUE);
  }
  if(any(is.na(unlist(tt)))){
    cat('Error: tt cannont contain NA/NaN entries\n')
    return(TRUE);
  } 
  ni_y = unlist(lapply(y,function(x) sum(!is.na(x))))
  if(all(ni_y == 1)){  
    cat("Error:FPCA is aborted because the data do not contain repeated measurements in y!\n"); 
    return(TRUE);    
  }
  ni_tt = unlist(lapply(tt,function(x) sum(!is.na(x))))
  if(all(ni_tt == 1)){  
    cat("Error:FPCA is aborted because the data do not contain repeated measurements in tt!\n"); 
    return(TRUE);    
  }   
  return(FALSE);
}

