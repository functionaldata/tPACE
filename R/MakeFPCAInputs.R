#' Format FPCA input
#'
#' Turn vector inputs to the list so they can be used in FPCA 
#' 
#' @param IDs  np-by-1 vector of subject IDs (Default: NULL)
#' @param tVec Either an np-by-1 vector of measurement times, or a p-by-1 vector corresponding to the common time support
#' @param yVec n-by-1 vector of measurements from the variable of interest, or a n-by-p matrix with each row corresponding to the dense observations.
#' @param na.rm logical indicating if NA should be omitted (Default: FALSE) 
#' @param sort logical indicating if Lt (and the correspoding Ly values) should be ensured to be sorted (Default: FALSE)
#' @param deduplicate logical indicating if the Lt should be ensured not to have within-subject duplicated values; the Ly values of repeated Lt values are averaged (Default: FALSE)
#' @return L list containing 3 lists each of length 'm', 'm' being the number of unique subject IDs
#' @export

MakeFPCAInputs <- function(IDs=NULL, tVec, yVec, na.rm=FALSE, sort=FALSE, deduplicate=FALSE){

  if ((!is.null(IDs) && any(is.na(IDs)) || any(is.na(tVec))) && 
      na.rm == FALSE) {
    stop('NAs exist in the IDs or tVec. Use na.rm=TRUE')
  }

  if( !is.null(IDs) ){ 
    if (na.rm) {
      dat <- na.omit(data.frame(IDs, tVec, yVec))
      IDs <- dat[, 'IDs']
      tVec <- dat[, 'tVec']
      yVec <- dat[, 'yVec']
    }
    uniqueIDs <- unique(IDs) 
    # Lt <- lapply( uniqueIDs, function(x) tVec[ which(IDs == x)])
    # Ly <- lapply( uniqueIDs, function(x) yVec[ which(IDs == x)])
    Lid <- as.list(uniqueIDs)
    Lt <- split(tVec, IDs, drop=TRUE)
    Lt <- Lt[match(as.character(uniqueIDs), names(Lt))]
    Ly <- split(yVec, IDs, drop=TRUE)
    Ly <- Ly[match(as.character(uniqueIDs), names(Ly))]
  } else if ( is.matrix(yVec) && is.null(IDs) && is.vector(tVec) ){
    if (ncol(yVec) != length(tVec)) {
      stop('columns of yVec does not correspond to tVec.')
    }
    Ly <- lapply( seq_len(nrow(yVec)), function(i) yVec[i,])
    Lt <- rep( list(tVec), dim(yVec)[1] )
    Lid <- as.list( 1:dim(yVec)[1])
  }
  if(sort){
    Ly = mapply( FUN = function(u1,u2) { ifelse(is.unsorted(u1), return( u2[order(u1)]), return(u2) ) }, u1 = Lt, u2 = Ly, SIMPLIFY = FALSE)
    Lt = lapply( Lt, function(u) { ifelse(is.unsorted(u), return( sort(u) ), return(u) ) } )
  }

  if(deduplicate){
    ids_with_duplicates <- which(unlist(lapply(Lt, function(x) length(x) != length(unique(x))))) 
    
    deduplicator <- function(x,y){
      xs_to_make_unique = sort(unique(x[duplicated(x)]))
      y_in_duplicated_xs = sapply(xs_to_make_unique, function(my_x) mean(y[my_x ==x]))
      new_x = c(x[!(x %in% xs_to_make_unique)], xs_to_make_unique) # unsorted
      new_y = c(y[!(x %in% xs_to_make_unique)], y_in_duplicated_xs) # unsorted
      return(list(sort(new_x), new_y[order(new_x)]))
    }
    
    if( length(ids_with_duplicates) > 0){
      for (my_id in ids_with_duplicates){
        no_dup_result = deduplicator(Lt[[my_id]], Ly[[my_id]])
        Lt[[my_id]] = no_dup_result[[1]]
        Ly[[my_id]] = no_dup_result[[2]]
      }
    }
  }
  
  L <- list( Lid = Lid, Ly = Ly, Lt = Lt)
  return(L)

}

