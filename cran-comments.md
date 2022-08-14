
## Test environments
* local R installation, macOS R 4.2.1
* local ubuntu 20.04, R 4.2.1
* win-builder (devel and release)

## R CMD check results
on macOS:
0 errors | 0 warnings | 0 notes

on windows:
0 errors | 0 warnings | 1 notes

Note indicates change in maintainer from Alvaro Gajardo to Yidong Zhou

on ubuntu one note is generated regarding installed package size of 1Mb or more in subfolder libs.

## Downstream dependencies
I have run R CMD check on ubuntu on downstream dependencies of fdapace: fdadensity,fdapaceShiny, fdaPOIFD, fgm, frechet, ftsa, KFPCA, LCox, mistat, SLFPCA and WRI.

All packages passed except LCox that produces an error on an example which does not seem to be related to a change in fdapace.

I have read and agree to all CRAN policies.
