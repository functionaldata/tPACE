IsRegular = function(t){
  
  # Check the if we have dense (2), or  dataType data with missing values (1) or sparse (0) data
  # t : n-by-1 list of vectors 
  
  tt = unlist(t);
  f = length(tt)/length(unique(tt))/length(t);
  if (f == 1){
    return('Dense');
  } else if(f > 0.75){
    return('RegularWithMV');
  } else {
    return('Sparse');
  }
}
