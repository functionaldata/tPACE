library(pracma)
source('../GetRawCov.R')
source('../mapX1d.R')
source('../IsRegular.R')
load('../data/dataForGetRawCov.RData')
GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',TRUE)  #Matches ML output
GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) #Matches ML output

y2 = list(1:10, 2:11)
t2 = list( 1:10, 1:10)

GetRawCov(y2,t2, sort(unique(unlist(t2))), seq(1.5,10.5, length.out=10) ,'Dense',TRUE) #Matches ML output
GetRawCov(y2,t2, sort(unique(unlist(t2))), seq(1.5,10.5, length.out=10) ,'Dense',FALSE) #Matches ML output