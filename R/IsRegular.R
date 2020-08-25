IsRegular = function(t){
  
  # Check the data type in terms of dense-sparse.
  # t : n-by-1 list of vectors 
  
  tt = unlist(t);
  f = length(tt)/length(unique(tt))/length(t);
  if (f == 1){
    if(length(unique(tt))<8){ #In case of low number of observations per subject
      return('Sparse');
    }
    else{
      return('Dense'); # for either regular and irregular data observed over the same time grid
    }
  }
  else{
    return('Sparse');
  }
}
