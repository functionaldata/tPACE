
## Test environments
* local R installation, macOS R 4.1.2
* local ubuntu 20.04, R 4.1.2
* win-builder (devel and release)

## R CMD check results
on macOS:
0 errors | 0 warnings | 0 notes

on ubuntu one note is generated regarding installed package size of 1Mb or more in subfolder libs. No notes were generated on windows.

## Downstream dependencies
I have run R CMD check on ubuntu on downstream dependencies of fdapace: fdadensity,fdapaceShiny, fdaPOIFD, fgm, frechet, ftsa, KFPCA, LCox, mistat, SLFPCA and WRI.

All packages passed.

I have read and agree to all CRAN policies.
