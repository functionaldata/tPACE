
library(devtools)
devtools::load_all("../../")
devtools::document("../../")
data=list.files(path="./",pattern="*.R")
data=data[-(which(data=="test_CreateStringingPlot.R"))]
data=data[-(which(data=="test_all.R"))]
for(i in 1:length(data)) {
  cat("Running test file: ", data[i], "\n")
  source(paste0("./",data[i]))
}
