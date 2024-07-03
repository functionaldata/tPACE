## Resubmission
This is a resubmission. In this version I have:
* Removed the old FLM function
* Updated FLM1 function and its documentation, which is now the FLM function
* Added FLMCI function to construct confidence intervals for functional linear models
* Changed @docType package to _PACKAGE in R/pkgname.R
* Added title and legend to CreateModeOfVarPlot function
* Updated figure titles in several plot functions (capitalize only the first word)


## Test environments
* local R installation, macOS R 4.4.1
* win-builder (devel and release)

## R CMD check results
on macOS:
0 errors | 0 warnings | 0 notes

on windows:
0 errors | 0 warnings | 0 notes

## Downstream dependencies
I have run R CMD check on macOS on downstream dependencies of fdapace: fdaconcur, fdadensity,fdapaceShiny, fdaPOIFD, fdarep, fgm, frechet, ftsa, KFPCA, longke, longsurr, MJMbamlss, mrct, SLFPCA and WRI.

All packages passed except fdaconcur, which produced an error due to the replacement of FLM1 with FLM. I have informed the maintainer of fdaconcur (we are in the same research group) about this change, and he will update the package to call `fdapace::FLM` in the next version.

I have read and agree to all CRAN policies.
