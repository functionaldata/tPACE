# # Based on https://github.com/hadley/testthat#integration-with-r-cmd-check

# Splitting test files into multiple ones so that each one runs within 10 mins,
# which is the limit on travis CI.
# This file contains the slow running tests.

library(testthat)
library(fdapace)
test_check("fdapace", filter='FClust|FSVD|FPCA|FVPA|funcReg', perl=TRUE)
