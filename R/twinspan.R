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
            ibegin = integer(max(mm, n)), idat = integer(ndat),
            PACKAGE = "twinspan")
    ibegin <- Z$ibegin
    idat <- Z$idat
    ## inflag: zero for species omitted as potential indicators
    inflag <- seq_len(n)
    if (!missing(noind))
        inflag[noind] <- 0
    ## Pseudospecies
    cutlevels <- as.integer(1000 * cutlevels + 0.5)
    nmax <- max(nlev * max(n, mm), 3 * (2^(levmax+1)-1))
    ## R cannot pass character vectors to Fortran, but we pass integer
    ## indices of names
    jname1 <- jname2 <- integer(nmax)
    jname1[seq_len(n)] <- jname2[seq_len(n)] <- seq_len(n)
    jnflag <- integer(nmax)
    jnflag[seq_len(n)] <- seq_len(n)
    iname1 <- iname2 <- integer(nmax)
    iname1[seq_len(mm)] <- iname2[seq_len(mm)] <- seq_len(mm)
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
    mmax <- max(mm, n)
    mmax <- max(mmax, nmax)
    jnam <- Z$jnam
    rrwt <- rep(1.0, mmax)
    ccwt <- rep(1.0, nmax)
    ccwt[jnam] <- lwgt[jnam] + TINY
    ## noind cases handled by inflag???
    indord <- rep(1L, nn)
    indpot <- Z$indpot
    jnflag <- Z$jnflag
    ## Call CLASS
    maxsam <- ndat # ??
    inddim <- indmax * (2^levmax-1) # dims for indicators
    Z <- .Fortran("class", mm=as.integer(mm), nn=as.integer(nn),
                  ndat=as.integer(ndat), mind=as.integer(indmax),
                  mmz=as.integer(MMZ), mms=as.integer(MMS),
                  ix=integer(mmax), iclass=integer(mmax),
                  iirow=integer(mmax), iaddr=ibegin,
                  indpot=indpot, indord=indord,
                  izone=integer(mmax), iy=Z$iy, jjcol=integer(nmax),
                  idat=Z$idat, indsig=integer(MMS),
                  ipict=integer(MMZ*MMS), x=double(nmax),
                  xx=double(nmax), rtot=double(mmax),
                  rrwt=as.double(rrwt), rowwgt=double(mmax),
                  y=double(nmax), yy=double(nmax),
                  ctot=double(nmax),
                  ccwt=as.double(ccwt), colwgt=double(nmax),
                  iname1=as.integer(iname1), iname2=as.integer(iname2),
                  jname1=Z$jname1, jname2=Z$jname2, jnam=Z$jnam,
                  x3=double(mmax), x4=double(mmax), x5=double(mmax),
                  lind=as.integer(lind), inflag=as.integer(inflag),
                  inlevmax=as.integer(levmax),
                  inmmin=as.integer(groupmin), eig = double(2^levmax-1),
                  indics = integer(inddim), isec = 1L, PACKAGE="twinspan")
    indics <- Z$indics
    dim(indics) <- c(indmax, 2^levmax-1)
    quadrat <- list(iclass = Z$iclass[seq_len(mm)], eig = Z$eig,
                    indicators = indics)
    ## species classification
    Y <- .Fortran("makejdat", mm=as.integer(mm), nn=as.integer(Z$nn),
                  nspec=as.integer(n), ndat=as.integer(Z$ndat),
                  nmax=as.integer(nmax), iaddr=as.integer(Z$iaddr),
                  idat=as.integer(Z$idat), indpot=as.integer(Z$indpot),
                  iclass=as.integer(Z$iclass), y=as.double(Z$y),
                  ccwt=as.double(Z$ccwt), rrwt=as.double(Z$rrwt),
                  jdat=integer(ndat),
                  PACKAGE="twinspan")
    Z <- .Fortran("class", nspec=as.integer(Y$nspec),
                  icmax=as.integer(3*max(Y$iclass)), ndat=as.integer(Y$ndat),
                  0L, mmz=as.integer(MMZ), mms=as.integer(MMS),  ix=Z$ix,
                  jnam=Z$jnam, iirow=Z$iirow, iaddr=Y$iaddr, indpot=Y$indpot,
                  indord=Z$indord, izone=Z$izone, iy=Z$iy, jjcol=Z$jjcol,
                  jdat=Y$jdat, indsig=Z$indsig, ipict=Z$ipict,
                  x=Z$x, xx=Z$xx, rtot=Z$rtot, rrwt=Y$rrwt, rowwgt=Z$rowwgt,
                  y=Y$y, yy=Z$yy, ctot=Z$ctot, ccwt=Y$ccwt, colwgt=Z$colwgt,
                  jname1=Z$jname1, jname2=Z$jname2, iname1=Z$iname1,
                  iname2=Z$iname2, inflag=Z$inflag, x3=Z$x3, x4=Z$x4, x5=Z$x5,
                  lind=Z$lind, jnflag=as.integer(jnflag),
                  inlevmax = Z$inlevmax, inmmin=Z$inmmin,
                  eig=double(2^levmax-1), indics = integer(inddim),
                  isec=2L,
                  PACKAGE="twinspan")
    species <- list(iclass = Z$jnam[seq_len(n)], eig = Z$eig)
    ## ordered index for quadrats and species
    quadrat$index <-
        .Fortran("clord", as.integer(mm), as.integer(levmax),
                       as.integer(quadrat$iclass), ix = integer(mm),
                       PACKAGE = "twinspan")$ix
    sindex <- .Fortran("clord", as.integer(n), as.integer(levmax),
                       as.integer(species$iclass), ix = integer(n),
                       PACKAGE = "twinspan")$ix
    irev <- 0
    rev <- .Fortran("revspec", nspec=as.integer(n), mm=as.integer(mm),
                    ndat=Z$ndat, ix=as.integer(quadrat$index),
                    iy=as.integer(sindex),
                    x=Z$x, y=Z$y, idat=Y$idat,
                    irev=as.integer(irev), PACKAGE="twinspan")$irev
    if (rev)
        sindex <- rev(sindex)
    species$index <- sindex
    ## out
    out <- list(call = match.call(), quadrat = quadrat, species = species)
    class(out) <- "twinspan"
    out
}
