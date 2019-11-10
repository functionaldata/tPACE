#' Number of eggs laid daily from medflies
#'
#' A dataset containing the eggs laid from 789 medflies (Mediterranean fruit flies,
#'  Ceratitis capitata) during the first 25 days of their lives. This is a subset of 
#' the dataset used by Carey at al. (1998); only flies that lived at least 25 days 
#' are included, i.e, at the end of the recording period [0,25] all flies are still alive.
#'
#' @name medfly25
#' @docType data
#' @format A data frame with 19725 rows and 3 variables:
#' \describe{
#' \item{ID}{: Medfly ID according to the original dataset}
#' \item{Days}{: Day of measurement}
#' \item{nEggs}{: Number of eggs laid at that particular day} 
#' \item{nEggsRemain}{: Remaining total number of eggs laid}
#' } 
#' @source \url{http://anson.ucdavis.edu/~mueller/data/medfly1000.html}
#' @references
#' {Carey, J.R., Liedo, P., Müller, H.G., Wang, J.L., Chiou, J.M. (1998). Relationship of age patterns of fecundity to mortality, longevity, and lifetime reproduction in a large cohort of Mediterranean fruit fly females. J. of Gerontology --Biological Sciences 53, 245-251. }
NULL
