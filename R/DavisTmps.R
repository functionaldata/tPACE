#' Monthly Air Temperature data in Davis, CA, USA as recorded in station GHCND:USC00042294
#'
#' A datasand containing the monthly air temperature data for Davis, CA, USA as this has been 
#' recorderd in NOAA's COOP station GHCND:USC00042294 from January 1893 to December 2015. All
#' temperature are reportedd in tenths of degrees Celsius and are derived by the daily maximum 
#' and minimum temperature during that month. Values of 9999 indicate missing readings.
#
#' @name DavisTemps
#' @docType data
#' @format A data frame with 1365 rows and 7 variables:
#' \describe{
#'   \item{EMXT}{: Extreme minimum temperature}
#'   \item{EMNT}{: Extreme maximum temperature}
#'   \item{MMNT}{: Monthly mean minimum temperature}
#'   \item{MMXT}{: Monthly mean maximum temperature}
#'   \item{MNTM}{: Monthly mean temperature}
#'   \item{Year}{: Year}
#'   \item{Month}{: Month; January:1, February:2, etc.}
#' } 
#' @source \url{http://www.ncdc.noaa.gov/cdo-web/datasets#GHCNDMS}
#' @references{Menne, M.J., I. Durre, R.S. Vose, B.E. Gleason, and T.G. Houston, 2012:  An overview of the Global Historical Climatology Network-Daily Database.  Journal of Atmospheric and Oceanic Technology, 29, 897-910, doi:10.1175/JTECH-D-11-00103.1.}
NULL
