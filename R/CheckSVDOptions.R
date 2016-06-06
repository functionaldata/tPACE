### CheckSVDOptions

CheckSVDOptions <- function(Ly1, Lt1, Ly2, Lt2, SVDoptns){
  
  if( (SVDoptns[['dataType1']]=='Sparse' && is.null(SVDoptns[['userMu1']])) ||
        (SVDoptns[['dataType2']]=='Sparse' && is.null(SVDoptns[['userMu2']])) ){
    stop('User specified mean function required for sparse functional data for cross covariance estimation.')
  }
  if(is.numeric(SVDoptns$methodSelectK)){
    if(SVDoptns$methodSelectK != round(SVDoptns$methodSelectK) || 
         SVDoptns$methodSelectK <= 0){
      stop('FSVD is aborted because the argument: methodSelectK is invalid!\n')
    }
  }

}