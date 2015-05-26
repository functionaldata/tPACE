# Translation of convertMuPhi.m
# One may rather just call mapX1d than using convertMuPhi.

convertMuPhi <- function(t, out1, mu, phi, regular) {

    if (regular == 'Dense') {
        muSub <- mapX1d(out1, mu, t[[1]])
        phiSub <- mapX1d(out1, phi, t[[1]])
    } else {
        muSub <- lapply(t, function(tt) mapX1d(out1, mu, tt))
        phiSub <- lapply(t, function(tt) mapX1d(out1, phi, tt))
    }
    
    return(list(muSub = muSub, phiSub = phiSub))
}