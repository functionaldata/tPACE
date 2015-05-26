
# Map (x,y) to (newx,newy)
# x, y : vectors of 1 * n
# newx : vector of 1 * m
# newy : vector of 1*m
# mapX1d <- function(x,y,newx){
  
  # if( is.vector(y) ){
    # return(y[is.element(x,newx)])
  # }else if(is.matrix(y)){
    # return(y[is.element(x,newx),])
  # }else{
    # warning('y cannot be empty!\n')
    # return(NaN)
  # }   
# }

mapX1d <- function(x, y, newx) {
    if (!all(newx %in% x)) 
        warning('Interpolation occured: you might want to increase the out1 coverage')
        
    if (min(newx) + 100 * .Machine$double.eps < min(x) || max(newx) > max(x) + 100 * .Machine$double.eps)
        warning('Extrapolation occured')

    newy <- approxExtrap(x, y, newx, method='linear')$y
    
    return(newy)
}