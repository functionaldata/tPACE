# # Based on https://github.com/hadley/testthat#integration-with-r-cmd-check

# Splitting test files into multiple ones so that each one runs within 10 mins,
# which is the limit on travis CI.
# This file contains the slow running tests.

library(testthat)
library(fdapace)

if (Sys.getenv('TRAVIS') != 'true') { 
test_check("fdapace", filter='FClust', perl=TRUE) #
test_check("fdapace", filter='FSVD', perl=TRUE) # over 10 min
test_check("fdapace", filter='FPCA', perl=TRUE) # over 10 min
test_check("fdapace", filter='FVPA', perl=TRUE) # over 10 min
test_check("fdapace", filter='FCReg', perl=TRUE) # over 10 min
test_check("fdapace", filter='FOptDes', perl=TRUE) #
test_check("fdapace", filter='GetCrCovYX', perl=TRUE) # over 10 min
test_check("fdapace", filter='selectK', perl=TRUE) #
test_check("fdapace", filter='WFDA', perl=TRUE) #
}
