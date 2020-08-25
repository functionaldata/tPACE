
library(devtools)
devtools::load_all("../../")
devtools::document("../../")
data=list.files(path="./",pattern="*.R")
data=data[-(which(data=="test_CreateStringingPlot.R"))]
data=data[-(which(data=="test_all.R"))]
sapply(1:length(data),function(i){source(paste0("./",data[i]))})