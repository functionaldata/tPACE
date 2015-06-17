devtools::load_all()
library(testthat)
options(error=recover)

trueLam <- 4 / ((2 * (1:50) - 1 ) * pi) ^ 2

set.seed(1)
n <- 100
pts <- seq(0, 1, by=0.05)
samp3 <- wiener(n, pts) + rnorm(n * length(pts), sd=0.1)
samp3 <- sparsify(samp3, pts, 10)
res <- FPCA(samp3$yList, samp3$tList, list(dataType='Sparse', useBins=TRUE))
res$lambda / trueLam[1:length(res$lambda)]
res$sigma2

createCorrPlot(res, 'Smoothed', FALSE)