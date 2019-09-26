# # Based on https://github.com/hadley/testthat#integration-with-r-cmd-check
# splitting test files into multiple ones so that each one runs within 10 mins,
# which is the limit on travis CI.
# This file contains the "fast" checks.

library(testthat)
library(fdapace)
test_check("fdapace", filter='_(?!FPCA)(?!FCReg)(?!FOptDes)(?!FClust)(?!FSVD)(?!FVPA)(?!GetCrCovYX)(?!selectK)(?!MultiFAM)(?!VCAM)(?!WFDA)', perl=TRUE)
# test_check("fdapace")

