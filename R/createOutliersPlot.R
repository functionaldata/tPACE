createOutliersPlot <- function(fpcaObj){

  if(length(fpcaObj$lambda) <2 ){
    stop("This plotting utility function needs at least two eigenfunctions")
    return(NULL)
  }
  
  xedge = 1.1 * max( abs(fpcaObj$xiEst[,1]))
  yedge = 1.1 * max( abs(fpcaObj$xiEst[,2]))  

  xlabelString = paste('FPC1 scores ', round(fpcaObj$FVE[1]), '%', sep=''   )
  ylabelString = paste('FPC2 scores ', round( diff( fpcaObj$FVE[1:2])), '%', sep=''   )

  plot(fpcaObj$xiEst[,1], fpcaObj$xiEst[,2], cex= .33, xlim = c(-xedge, xedge), ylim =c(-yedge, yedge), 
       pch=10, xlab= xlabelString, ylab=ylabelString )
  abline(v=(seq(-xedge, xedge, length.out=  21)), col="lightgray", lty="dotted")
  abline(h=(seq(-yedge, yedge, length.out=  21)), col="lightgray", lty="dotted")

  varXi1Xi2 = var(cbind(fpcaObj$xiEst[,1], fpcaObj$xiEst[,2]))
  muXi1Xi2 =  apply(cbind(fpcaObj$xiEst[,1], fpcaObj$xiEst[,2]),2,mean)

  lines(ellipse::ellipse( varXi1Xi2, centre= muXi1Xi2, level=.5), col=2)
  lines(ellipse::ellipse( varXi1Xi2, centre= muXi1Xi2, level=.95), col=3)
  lines(ellipse::ellipse( varXi1Xi2, centre= muXi1Xi2, level=.99), col=4)
  lines(ellipse::ellipse( varXi1Xi2, centre= muXi1Xi2, level=.999), col=5)   
  legend(legend= c('0.500', '0.950', '0.990', '0.999'), x='topright', col=2:5, lwd=2)

}
