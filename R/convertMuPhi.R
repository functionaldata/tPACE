# Translation of convertMuPhi.m
# One may rather just call mapX1d than using convertMuPhi.

convertMuPhi <- function(t, obsGrid, mu, phi, dataType) {

    if (dataType == 'Dense') {
        muSub <- mapX1d(obsGrid, mu, t[[1]])
        phiSub <- mapX1d(obsGrid, phi, t[[1]])
    } else {
        muSub <- lapply(t, function(tt) mapX1d(obsGrid, mu, tt))
        phiSub <- lapply(t, function(tt) mapX1d(obsGrid, phi, tt))
    }
    
    return(list(muSub = muSub, phiSub = phiSub))
}