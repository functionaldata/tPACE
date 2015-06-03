## To be merged into FPCA.R

CreateOptions = function(bwmu = NULL, bwmuGcv = NULL, bwuserCov = NULL, bwuserCovGcv = NULL, 
    ntest1 = NULL, ngrid1 = NULL, selectionMethod = NULL, FVEthreshold = NULL,
    maxK = NULL, dataType = NULL, error = NULL, nRegGrid = NULL,
    method = NULL, shrink = NULL, newdata = NULL, kernel = NULL, 
    numBins = NULL, yname = NULL, screePlot = NULL, designPlot = NULL, 
    corrPlot = NULL,     rho = NULL, verbose = NULL, userMu = NULL, userCov = NULL, methodCov = NULL,
    methodMu = NULL, outPercent = NULL, useBinnedData = NULL, rotationCut = NULL){ 

 return( list(bwmu = bwmu, bwmuGcv = bwmuGcv, bwuserCov = bwuserCov, bwuserCovGcv = bwuserCovGcv,
          ntest1 = ntest1, ngrid1 = ngrid1, selectionMethod = selectionMethod, FVEthreshold = FVEthreshold,
          maxK = maxK, dataType = dataType, error = error, nRegGrid = nRegGrid, 
          method = method, shrink = shrink, newdata = newdata, kernel = kernel, corrPlot = corrPlot,	
          numBins = numBins, yname = yname, screePlot = screePlot, designPlot = designPlot, rho = rho, rotationCut = rotationCut,
          verbose = verbose, userMu = userMu, userCov = userCov, methodMu= methodMu,  methodCov= methodCov, outPercent = outPercent, useBinnedData = useBinnedData, rotationCut = rotationCut) )
}
