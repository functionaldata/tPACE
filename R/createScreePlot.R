
createScreePlot <-function(ys){ 
 
  dfbar <-barplot( rep(NA,length(ys)), ylim=c(0,105) ,axes=FALSE, xlab='Number of components', ylab='Fraction of Variance Explained', names.arg=1:length(ys), main="Screeplot")


  abline(h=(seq(0,100,5)), col="lightgray", lty="dotted")
  barplot(c(ys[1], diff(ys)), add = TRUE )
  lines(dfbar, y= ys, col='red')
  points(dfbar, y= ys, col='red')
  legend("right", "Cummul. FVE", col='red', lty=1, pch=1, bty='n') 
   
}

