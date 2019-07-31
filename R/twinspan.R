### R function to call needed parts of TWINSPAN Fortran functions. The
### original TWINSPAN was split into separate files for each
### subroutine, and this R function will call those subroutines as
### needed, but R code is used instead of FORTRAN main. This allows us
### to make the TWINSPAN work like an ordinary R function.

#' Two-Way Indicator Species Analysis
#'
#' Two-Way Indicator Analysis (TWINSPAN) is a divisive classification
#' method that works by splitting first Correspondence-Analysis into
#' two, and then recursively working with each of this split
#' subsets. The function is actually much more complicated: see
#' Details.
#'
#' @param x Input data, usually a species community data set where
#'     columns give the species and rows the sampling units.
#'
#' @useDynLib twinspan
#'
#' @export
`twinspan` <-
    function(x)
{
    x <- as.matrix(x)
    n <- ncol(x) # no. of species
    mm <- nrow(x) # no. of SUs
    nid <- sum(x > 0) # number of positive items
    ## translate data to the internal sparse format
    Z <- .C("data2hill", as.double(x), mm = as.integer(mm),
            n = as.integer(n), nid = as.integer(nid),
            ibegin = integer(mm), iend = integer(mm),
            idat = integer(nid), qidat = double(nid),
            PACKAGE = "twinspan")
    ibegin <- Z$ibegin
    iend <- Z$iend
    idat <- Z$idat
    qdat <- Z$qidat
    Z
}
