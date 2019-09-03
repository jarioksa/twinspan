### R function to call needed parts of TWINSPAN Fortran functions. The
### original TWINSPAN was split into separate files for each
### subroutine, and this R function will call those subroutines as
### needed, but R code is used instead of FORTRAN main. This allows us
### to make the TWINSPAN work like an ordinary R function.

#' Two-Way Indicator Species Analysis
#'
#' Two-Way Indicator Analysis (TWINSPAN) is a divisive classification
#' method that works by splitting first Correspondence Analysis into
#' two classes, and then recursively working with each split
#' subset. The function is actually much more complicated: see
#' Details.
#'
#' @param x Input data, usually a species community data set where
#'     columns give the species and rows the sampling units.
#' @param cutlevels Cut levels used to split quantitative data into
#'     binary pseudospecies.
#' @param indmax Maximum number of indicators for division.
#' @param groupmin Minimum group size for division.
#' @param lind Weights for levels of pseudospecies. For example
#'     indicator potentials \code{c(1, 0, 0,1, 0)} signify that
#'     pseudospecies at levels 1 and 4 can be used as indicators, but
#'     that those at other levels cannot. In the default case, all
#'     species are available.
#' @param levmax Maximum level of divisions.
#' @param lwgt Weights for the levels of pseudospecies. For example
#'     weights \code{c(1, 2, 2, 2)} signify that pseudospecies
#'     corresponding to 3 higher cut levels are to be given twice the
#'     weight of pseudospecies at the lowest level.
#' @param noind Numbers (indices) of species that you wish to omit
#'     from list of potential indicators. Species omitted from this
#'     list are used in the calculation, but cannot appear as
#'     indicators.
#'
#' @useDynLib twinspan
#'
#' @export
`twinspan` <-
    function(x, cutlevels = c(0,2,5,10,20), indmax = 7,
             groupmin = 5, levmax = 6,
             lind, lwgt, noind)
{
    ## handle arguments
    if (missing(lind))
        lind <- rep(1, length(cutlevels))
    if (missing(lwgt))
        lwgt <- rep(1, length(cutlevels))
    ipunch = 0L ## never write to a file
    ## set internal constants and their derived quantities
    MZCRIT <- 8L
    MZOUT <- 4L
    MZ <- as.integer(MZCRIT + 2*MZOUT)
    MMZ <- MZ + 4L
    MMS <- as.integer(indmax + 4L)
    TINY <- 1e-5
    ##
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
    ## inflag: zero for species omitted as potential indicators
    inflag <- seq_len(n)
    if (!missing(noind))
        inflag[noind] <- 0
    ## Pseudospecies
    cutlevels <- as.integer(1000 * cutlevels + 0.5)
    nmax <- nlev * max(n, mm)
    ## R cannot pass character vectors to Fortran, but we pass integer
    ## indices of names
    jname1 <- jname2 <- integer(nmax)
    jname1[seq_len(n)] <- jname2[seq_len(n)] <- seq_len(n)
    jnflag <- integer(nmax)
    jnflag[seq_len(n)] <- seq_len(n)
    iname1 <- iname2 <- seq_len(mm)
    Z <- .Fortran("pseudo", mm = as.integer(mm), nn = as.integer(n),
                  nmax = as.integer(nmax), nl = as.integer(nlev),
                  ndat = as.integer(ndat), nspec = as.integer(nmax),
                  idat = as.integer(idat), lcut = as.integer(cutlevels),
                  jnflag = as.integer(jnflag),
                  jname1 = as.integer(jname1), jname2 = as.integer(jname2),
                  jnam = integer(nmax), indpot = integer(nmax),
                  iy = integer(nmax), PACKAGE = "twinspan")
    nn <- Z$nn
    mm <- Z$mm
    jnam <- Z$jnam
    rrwt <- rep(1.0, nmax)
    ccwt <- rep(1.0, nmax)
    ccwt[jnam] <- lwgt[jnam] + TINY
    ## noind cases handled by inflag???
    indord <- rep(1L, nn)
    indpot <- Z$indpot
    ## Call CLASS
    maxsam <- ndat # ??
    Z <- .Fortran("class", mm=as.integer(mm), nn=as.integer(nn),
                  ndat=as.integer(ndat), mind=as.integer(indmax),
                  mmz=as.integer(MMZ), mms=as.integer(MMS),
                  ix=integer(mm), iclass=integer(mm),
                  iirow=integer(mm), iaddr=ibegin,
                  indpot=indpot, indord=indord,
                  izone=integer(mm), iy=Z$iy, jjcol=integer(nmax),
                  idat=Z$idat, indsig=integer(indmax),
                  ipict=integer(104*25), x=double(mm),
                  xx=double(mm), rtot=double(mm),
                  rrwt=as.double(rrwt), rowwgt=double(mm),
                  y=double(nn), yy=double(nn), ctot=double(nn),
                  ccwt=as.double(ccwt), colwgt=double(nn),
                  iname1=as.integer(iname1), iname2=as.integer(iname2),
                  jname1=Z$jname1, jname2=Z$jname2, jnam=Z$jnam,
                  x3=double(mm), x4=double(mm), x5=double(mm),
                  lind=as.integer(lind), inflag=as.integer(inflag),
                  inlevmax=as.integer(levmax),
                  inmmin=as.integer(groupmin),
                  PACKAGE="twinspan")
    ## species classification
    Y <- .Fortran("makejdat", mm=as.integer(mm), nn=as.integer(Z$nn),
                  nspec=as.integer(n), ndat=as.integer(Z$ndat),
                  nmax=as.integer(nmax), iaddr=as.integer(Z$iaddr),
                  idat=as.integer(Z$idat), indpot=as.integer(Z$indpot),
                  iclass=as.integer(Z$iclass), y=as.double(Z$y),
                  ccwt=as.double(Z$ccwt), rrwt=as.double(Z$rrwt),
                  jdat=integer(ndat),
                  PACKAGE="twinspan")

## out
    Z$call <- match.call()
    class(Z) <- "twinspan"
    Z
}
