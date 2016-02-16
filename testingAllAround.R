
rm(list=ls())
setwd('~/git_projects/tPACE/')
devtools::load_all()
A = list.files('tests/testthat/') 

for (i in 1:length(A)){
  print(A[i])
  eval(parse(text=paste(collapse = '',c("source('tests/testthat/",A[i] ,"')") )))

}


