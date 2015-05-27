# This version takes a 1 by bwlen candidate bandwidth vector input, rather than finding out the candidates itself.
gcvlwls2d1 <- function(t_all, ngrid, regular, error, kern, rcov,
                      verbose=FALSE, bwCandidates) {
# get the regression model
    r <- max(t_all) - min(t_all)

    gcvScores <- sapply(bwCandidates, function(bw) gcv(...))


    if(all((is.infinite(gcvScores)))){
        error("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")      
    }

    bInd = which(gcvScores == min(gcvScores));
    bScr = gcvScores[bInd][1]
    bOpt = max(bwCandidates[bInd]);

    if( bOpt == r) {
      warning("data is too sparse, optimal bandwidth includes all the data!You may want to change to Gaussian kernel!\n")
    }
    bOptList <- list( 'bOpt' = bOpt, 'bScore' = bScr) 
    return( bOptList)
    
}

# function [bw_xcov,gcv]=gcv_mullwlsn(t,ngrid,regular,error,kernel,rcov, verbose)
