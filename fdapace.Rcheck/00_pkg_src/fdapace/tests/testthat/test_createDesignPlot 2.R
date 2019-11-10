library(testthat)
# devtools::load_all()

# Uniform time points
## a speed test
set.seed(1)
n <- 1e3
sparsity <- 1:5
Lt <- replicate(n, runif(sample(sparsity, 1)), simplify=FALSE)
obsGrid <- sort(unique(unlist(Lt)))
system.time(
CreateDesignPlot(Lt, obsGrid, isColorPlot=TRUE)
)

# ... are passed in 
set.seed(1)
n <- 5e2
sparsity <- 1:5
Lt <- replicate(n, round(runif(sample(sparsity, 1)), 2), simplify=FALSE)
obsGrid <- sort(unique(unlist(Lt)))
CreateDesignPlot(Lt, obsGrid, isColorPlot=TRUE, pch=1, cex=1, xlab='XX', ylab='YY')
