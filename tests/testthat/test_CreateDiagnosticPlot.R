# devtools::load_all()

viewPdf <- FALSE
set.seed(1)
n <- 100
M <- 51
pts <- seq(0, 1, length.out=M)
mu <- rep(0, length(pts))
sampDense <- Wiener(n, pts)
samp <- Sparsify(sampDense, pts, M)
res <- FPCA(samp$Ly, samp$Lt, list(error=TRUE, FVEthreshold=1, dataType='Dense', plot=TRUE))

if (viewPdf) {
  pdf('tmp.pdf')
}

plot(res)
CreateDesignPlot(samp$Lt)
plot(res, addLegend=FALSE)
CreateDesignPlot(samp$Lt, addLegend=FALSE)

if (viewPdf) {
  dev.off()
  system('open tmp.pdf')
  file.remove('tmp.pdf')
}
