
# Map (x,y) to (newx,newy)
# x, y : vectors of 1 * n
# newx : vector of 1 * m
# newy : vector of 1*m
mapX1d <- function(x,y,newx){
  
  if( is.vector(y) ){
    return(y[is.element(x,newx)])
  }else if(is.matrix(y)){
    return(y[is.element(x,newx),])
  }else{
    warning('y cannot be empty!\n')
    return(NaN)
  }   
}
