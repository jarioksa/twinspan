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
#' @param cutlevels Cut levels used to split quantitative data into
#'     binary pseudospecies.
#'
#' @useDynLib twinspan
#'
#' @export
`twinspan` <-
    function(x, cutlevels = c(0,2,5,10,20))
{
    nlev <- length(cutlevels)
    x <- as.matrix(x)
    n <- ncol(x) # no. of species
    mm <- nrow(x) # no. of SUs
    nid <- 2 * sum(x > 0) + mm # length of idat
    ndat <- 2 * nlev * nid
    ## translate data to the internal sparse format
    Z <- .C("data2hill", as.double(x), mm = as.integer(mm),
            n = as.integer(n), nid = as.integer(nid),
            ibegin = integer(mm), idat = integer(ndat),
            PACKAGE = "twinspan")
    ibegin <- Z$ibegin
    idat <- Z$idat
    ## Pseudospecies
    cutlevels <- as.integer(1000 * cutlevels + 0.5)
    nmax <- nlev * n
    ## R cannot pass character vectors to Fortran, but we pass integer
    ## indices of names
    jname1 <- jname2 <- integer(nmax)
    jname1[seq_len(n)] <- jname2[seq_len(n)] <- seq_len(n)
    jnflag <- integer(nmax)
    jnflag[seq_len(n)] <- seq_len(n)
    Z <- .Fortran("pseudo", mm = as.integer(mm), nn = as.integer(n),
                  nmax = as.integer(nmax), nl = as.integer(nlev),
                  ndat = as.integer(ndat), nspec = as.integer(nmax),
                  idat = as.integer(idat), lcut = as.integer(cutlevels),
                  jnflag = as.integer(jnflag),
                  jname1 = as.integer(jname1), jname2 = as.integer(jname2),
                  jnam = integer(nmax), indpot = integer(nmax),
                  iy = integer(nmax), PACKAGE = "twinspan")
    Z
}
