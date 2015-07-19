CreateEigenfncts = function(tt, lint = c(0,1), numPhi = 2){
  # Former xeig  
  
  if(length(lint) == 1){
    lb = 0;
    Tc = lint;
  } else {
    lb = lint[1];
    Tc = lint[2];
  }
  
  id = which(tt < lb | tt > Tc)
  if(length(id) > 0){
    cat("tt must be in [", lb , ",", Tc, "]. Invalid elements in tt are removed!\n")
    tt = tt[-id]
  }
  
  phi = matrix(NA,numPhi, length(tt));
  
  if(numPhi == 1){
    phi = -sqrt(2/Tc)*cos(2*pi*tt/Tc);
    phi = matrix(phi, 1, length(phi));
    
  } else if(numPhi == 2){
    phi[1,] = -sqrt(2/Tc)*cos(2*pi*tt/Tc);
    phi[2,] = sqrt(2/Tc)*sin(2*pi*tt/Tc);    
  } else {
    phi = matrix(0, numPhi, length(tt));
    id = 1:numPhi;
    oddID = id%%2;
    oddFactor = 1:sum(oddID);
    evenID = oddID == 0;
    evenFactor = 1:sum(evenID);
    
    phiOdd = matrix(0, sum(oddID), length(tt));
    phiEven = matrix(0, sum(evenID), length(tt));
    
    for(i in 1:sum(oddID)){
      phiOdd[i,] = -sqrt(2/Tc)*cos(2*oddFactor[i]*pi*tt/Tc);
    }
    phi[which(oddID == 1),] = phiOdd
    
    for(i in 1:sum(evenID)){
      phiEven[i,] = sqrt(2/Tc)*sin(2*evenFactor[i]*pi*tt/Tc);
    }
    phi[which(evenID == 1),] = phiEven;
  }
  
  return(phi)   
  
}
