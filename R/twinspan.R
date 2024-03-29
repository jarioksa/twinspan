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
#' subset. The current function is based on and uses much of the
#' original FORTRAN code of the original TWINSPAN (Hill
#' 1979). \code{twinspan} is the main function of this package, but it
#' works silently and prints very little information: you must use
#' separate support functions to extract various aspects of the
#' result.
#'
#' \code{twinspan} may not print anything when it runs, but it will
#' return its result that you should save for later use. The functions
#' that reproduce most of the traditional printout are
#' \code{\link{summary.twinspan}} and \code{\link{twintable}}. The
#' \code{summary} prints the history of divisions with eigenvalues,
#' signed indicator pseudospecies and the threshold of indicator
#' scores for division, and for terminal groups it prints the group
#' size and group members (quadrats, species). Function
#' \code{twintable} prints the classified community table. In
#' addition, \code{plot.twinspan} shows the dendrogram corresponding
#' to the \code{summary} with division numbers and sizes and id number
#' of terminal groups. Function \code{\link{image.twinspan}} provides
#' a graphical overview of major structure of classification as a
#' prelude to \code{twintable}. With function
#' \code{\link{as.dendrogram.twinspan}} it is possible to construct a
#' \code{\link{dendrogram}} of complete classification down to final
#' units (quadrats, species), and \code{\link{as.hclust.twinspan}}
#' constructs an \code{\link{hclust}} tree down to final groups.
#'
#' The \code{twinspan} function performs the classic TWINSPAN with
#' fixed levels of hierarchy, but with other functions in this package
#' it is also possible to perform the modified method of Roleček et
#' al. (2009): the topography of the tree is the same, but the
#' division heights are based on the hetergeneity of the group, and
#' the groups are extracted in the order of their heterogeneity. The
#' measure of the heterogeneity is total standardized chi-square (or
#' inertia) of the divided group. This is based on the same matrix
#' as used internally in \code{twinspan} code (see support functions
#' \code{\link{twintotalchi}}, \code{\link{twin2stack}} and
#' \code{\link{twin2specstack}}).
#'
#' The classification at any level of division can be extracted with
#' \code{\link{cut.twinspan}}, and the most heterogenous groups with
#' \code{\link{cuth}}. Function \code{\link{predict.twinspan}}
#' provides a similar classification vector, but based on indicator
#' pseudospecies, and can be used also with new data that was not used
#' in \code{twinspan}. These two classifications are often in
#' conflict, and \code{\link{misclassified}} will detect those cases
#' and the divisions where the two classifications diverged. Function
#' \code{\link{eigenvals}} extracts the eigenvalues of divisions, and
#' \code{\link{twintotalchi}} finds the \dQuote{sum of all
#' eigenvalues} or the standardized chi-square of each division or
#' final group.
#'
#' Function \code{\link{twinsform}} transforms the data similarly as
#' \code{twinspan} and can be used to reproduce the results of any
#' single division. Functions \code{\link{twin2mat}} and
#' \code{\link{twin2stack}} extract the internal data matrices in
#' standard \R{} format from the \code{twinspan} result.
#'
#' @section Method:
#'
#' TWINSPAN is very complicated and has several obscure details, and
#' it will not be explained in details in this manual, but you should
#' consult the source code or literature sources. Hill (1979) is the
#' most authorative source, but may be difficult to find. Kent & Coker
#' (1991) do a great job in explaining the method, including many
#' obscure details.
#'
#' A strong simplification (but often sufficient to understand the
#' basic principles) is that TWINSPAN is a divisive clustering based
#' on splitting first correspondence analysis axis, and applying the
#' same method recursively for resulting classes. The same method is
#' used first for quadrats and then for species. In addition, it finds
#' the species abundance levels (called \sQuote{pseudospecies}) that
#' best indicate these divisions facilitating ecologist's
#' understanding of classes. Species classification is performed so
#' that it best corresponds to the previous quadrat classification
#' also with species composition.
#'
#' The following details are more technical. The analysis starts with
#' splitting species abundance data into discrete abundance levels
#' called pseudospecies. With these the function constructs a stacked
#' binary matrix with values of 0 and 1, where value 1 means that
#' species occurs at given threshold (called \sQuote{cut level}) in a
#' quadrat. Then the pseudospecies (abundance levels) that occur in
#' fewer than 1/5 of quadrats are downweighted so that their presences
#' (values 1) are reduced linearly towards minimimum value of 0.01
#' according to their frequencies. This will reduce their impact in
#' correspondence analysis which is regarded as being sensitive to
#' rare species. Each cut level of the species is downweighted
#' independently. The first axis of correspondence analysis is found
#' for the downweighted data. This initial step can be reproduced with
#' the help of functions \code{\link{twinsform}} or
#' \code{\link{twin2stack}}.  However, the division is not based on
#' this step only. Next the method finds the best indicator
#' pseudospecies for the division. Further, it polarizes the
#' ordination by using indicator scores for all species to find the
#' final classes for quadrats. It does not mechanically just
#' split the axis in the middle, but it finds the cutpoint so that
#' indicator scores from the indicator pseudospecies and final split
#' are as concordant as possible. Then the analysis is repeated for
#' both resulting groups, including downweighting within the subset of
#' quadrats.
#'
#' After quadrat classification, TWINSPAN constructs species data
#' which are completely different from the data used in quadrat
#' classification.  Species values depend on their ability to
#' discriminate quadrat classes at any level of classification. Then
#' species are classified in the same way and with the same code as
#' the quadrats. In this way species classification is concordant with
#' quadrat classification, and good indicators of quadrat classes are
#' grouped together. The species classification can be reproduced with
#' function \code{\link{twin2specstack}} which also provides a more
#' detailed description of the data structure used at this stage.
#'
#' Function \code{twinspan} performs only the classical TWINSPAN, but
#' with support functions the modified method of Roleček et al. (2009)
#' can be performed (see \code{\link{cuth}},
#' \code{\link{as.hclust.twinspan}}).
#'
#' @encoding UTF-8
#'
#' @references
#'
#' Hill, M.O. (1979). \emph{TWINSPAN - a FORTRAN program for arranging
#' multivariate data in an ordered two-way table by classification of
#' individuals and attributes.} Cornell Univ., Dept of Ecology and
#' Systematics.
#'
#' Kent, M. & Coker, P. (1992) \emph{Vegetation description and
#' analysis: A practical approach.} John Wiley & Sons.
#'
#' Roleček, J, Tichý, L., Zelený, D. & Chytrý, M. (2009). Modified
#' TWINSPAN classification in which the hierarchy respects cluster
#' heterogeneity. \emph{J Veg Sci} 20: 596--602.
#'
#' @return
#'
#' Function returns an object of class \code{"twinspan"}, with
#' following items: \describe{
#'
#' \item{call}{Function call.}
#'
#' \item{cutlevels}{Defined cutlevels. These will be used in
#'   \code{\link{predict.twinspan}}.}
#'
#' \item{levelmax}{Maximum level depth of divisions. The divisions
#' will end when this depth is achieved.}
#'
#' \item{nspecies}{Number of species.}
#'
#' \item{nquadrat}{Number of quadrats.}
#'
#' \item{idat}{Pseudospecies data in the internal format used in
#' \code{twinspan}.  Functions \code{\link{twin2mat}} and
#' \code{\link{twin2stack}} can change this into more usable format.}
#'
#' \item{quadrat}{Results for quadrats (described below).}
#' \item{species}{Results for species (described below).}
#' }
#'
#' The results of the analysis are stored in items \code{quadrat} and
#' \code{species} with similar structure, but \code{species} has only
#' items \code{iclass}, \code{eig}, \code{labels} and \code{index} of
#' the following:
#'
#' \describe{
#'
#' \item{iclass}{ID numbers of final classes at the lowest level of
#' hierarchy. These can be extracted with \code{\link{cut.twinspan}}
#' which also can transform these to any higher level of hierarchy.}
#'
#' \item{eig}{Eigenvalues of divisions.}
#'
#' \item{labels}{Name labels of units (species, quadrats).}
#'
#' \item{index}{An index to order the units in \code{twinspan}
#' displays, e.g., in \code{\link{twintable}}.}
#'
#' \item{indicators}{A matrix of dimensions maximum number of
#' indicators \eqn{\times}{times} maximum number of divisions giving the
#' signed indices of indicator species for the division. These are
#' shown in labels in \code{\link{summary.twinspan}} and used by
#' \code{\link{predict.twinspan}}.}
#'
#' \item{positivelimit}{Lowest value of indicator score for positive
#' group in a division.}
#'
#' \item{indlabels}{Labels of pseudospecies.}
#'
#' \item{pseudo2species}{Index from pseudospecies to the corresponding
#' species.}
#' }
#'
#' @examples
#'
#' data(ahti)
#' ## default cut levels
#' (tw <- twinspan(ahti))
#' ## visual look at the divisions and group numbers
#' plot(tw)
#' ## Braun-Blanquet scale
#' (twb <- twinspan(ahti, cutlevels = c(0, 0.1, 1, 5, 25, 50, 75)))
#' plot(twb)
#' ## compare confusion
#' table(cut(tw, level=3), cut(twb, level=3))
#' ## modified method of Roleček et al. (2009)
#' plot(twb, height="chi", main = "Rolecek tree")
#' ## compare against the default by hierarchy levels
#' table(cuth(twb, ngroups=8), cut(twb, level=3))
#'
#'
#' @param x Input data, usually a species community data set where
#'     columns give the species and rows the sampling units.
#' @param cutlevels Cut levels used to split quantitative data into
#'     binary pseudospecies. Max of 9 cutlevels can be used.
#' @param indmax Maximum number of indicators for division (15 or less).
#' @param groupmin Minimum group size for division (2 or larger).
#' @param lind Weights for levels of pseudospecies. For example
#'     indicator potentials \code{c(1, 0, 0,1, 0)} signify that
#'     pseudospecies at levels 1 and 4 can be used as indicators, but
#'     that those at other levels cannot. In the default case, all
#'     species are available.
#' @param levmax Maximum depth of levels of divisions (15 or less).
#' @param lwgt Weights for the levels of pseudospecies. For example
#'     weights \code{c(1, 2, 2, 2)} signify that pseudospecies
#'     corresponding to 3 higher cut levels are to be given twice the
#'     weight of pseudospecies at the lowest level.
#' @param noind Numbers (indices) of species that you wish to omit
#'     from list of potential indicators. Species omitted from this
#'     list are used in the calculation, but cannot appear as
#'     indicators.
#'
#' @useDynLib twinspan, .registration=TRUE
#'
#' @export
`twinspan` <-
    function(x, cutlevels = c(0,2,5,10,20), indmax = 7,
             groupmin = 5, levmax = 6,
             lind, lwgt, noind)
{
    ## handle arguments
    if (length(cutlevels) > 9)
        stop("max 9 cutlevels accepted")
    if (groupmin < 2)
        stop("groupmin must be at least 2")
    if (indmax > 15)
        stop("indmax must be 15 or less")
    if (levmax > 15)
        stop("levmax must be 15 or less")
    if (missing(lind))
        lind <- rep(1, length(cutlevels))
    else
        if (length(lind) != length(cutlevels))
            stop("lind must have same length as cutlevels")
    if (missing(lwgt))
        lwgt <- rep(1, length(cutlevels))
    else
        if (length(lwgt) != length(cutlevels))
            stop("lwgt must have same length as cutlevels")
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
    ## remove empty species (if any)
    csum <- colSums(x)
    if (any(csum == 0)) {
        warning("some empty species were removed")
        x <- x[, csum>0, drop=FALSE]
    }
    ## stop with empty quadrats: twinspan could manage them, but many
    ## support functions should be re-written: better have users to
    ## remove empty quadrats.
    if (any(rowSums(x) <= 0))
        stop("you have empty quadrats: please remove them")
    n <- ncol(x) # no. of species
    mm <- nrow(x) # no. of SUs
    nid <- 2 * sum(x > 0) + mm # length of idat
    ndat <- 2 * nlev * nid
    ## translate data to the internal sparse format
    Z <- .C("data2hill", as.double(x), mm = as.integer(mm),
            n = as.integer(n), nid = as.integer(nid),
            ibegin = integer(max(mm, n)), idat = integer(ndat))
    ibegin <- Z$ibegin
    idat <- Z$idat
    ## inflag: zero for species omitted as potential indicators
    inflag <- seq_len(n)
    if (!missing(noind))
        inflag[noind] <- 0
    ## Pseudospecies
    cut1000 <- as.integer(1000 * cutlevels + 0.5)
    nmax <- max(nlev * max(n, mm), 3 * (2^(levmax+1)-1))
    ## R cannot pass character vectors to Fortran, but we pass integer
    ## indices of names
    jname1 <- integer(nmax)
    jname1[seq_len(n)] <- seq_len(n)
    jnflag <- integer(nmax)
    jnflag[seq_len(n)] <- seq_len(n)
    Z <- .Fortran("pseudo", mm = as.integer(mm), nn = as.integer(n),
                  nmax = as.integer(nmax), nl = as.integer(nlev),
                  ndat = as.integer(ndat), nspec = as.integer(nmax),
                  idat = as.integer(idat), lcut = as.integer(cut1000),
                  jnflag = as.integer(jnflag),
                  jname1 = as.integer(jname1),
                  jnam = integer(nmax), indpot = integer(nmax),
                  iy = integer(nmax))
    ## data & names
    idat <- Z$idat[1:which(Z$idat==-1)[mm]]
    rnames <- rownames(x, do.NULL = FALSE, prefix = "q")
    cnames <- colnames(x, do.NULL = FALSE, prefix = "sp")
    k <- Z$jname1 > 0
    pnames <- paste0(cnames[Z$jname1[k]], Z$jnam[k])
    ## some other re-set arguments
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
                  jname1=Z$jname1, jnam=Z$jnam,
                  x3=double(mmax), x4=double(mmax), x5=double(mmax),
                  lind=as.integer(lind),
                  inlevmax=as.integer(levmax),
                  inmmin=as.integer(groupmin), eig = double(2^levmax-1),
                  indics = integer(inddim), limpos = integer(2^levmax-1),
                  isec = 1L)
    indics <- Z$indics
    dim(indics) <- c(indmax, 2^levmax-1)
    quadrat <- list(iclass = Z$iclass[seq_len(mm)], eig = Z$eig,
                    indicators = indics, positivelimit = Z$limpos,
                    labels = rnames, indlabels = pnames,
                    pseudo2species = jnflag[jnflag > 0])
    ## species classification: build a matrix of preferences for
    ## terminal clusters & internal divisions: see twin2specstack()
    iclass <- quadrat$iclass
    spmat <- matrix(0, max(iclass), n)
    ## species totals in each terminal group
    icl <- iclass[1]
    k <- 1
    spsum <- numeric(n)
    for(i in seq_along(idat)) {
        ## new quadrat?
        if (idat[i] == -1) {
            k <- k + 1
            icl <- iclass[k]
            next
        }
        ## new species?
        if (idat[i] <= n) {
            j <- idat[i]
            spsum[j] <- spsum[j] + 1L
        }
        spmat[icl, j] <- spmat[icl, j] + 1L
    }
    ## accumulate species totals for internal divisions
    for (i in max(iclass):2) {
        spmat[i %/% 2,] <- spmat[i %/% 2,] + spmat[i,]
    }
    ## class totals
    clsum <- rowSums(spmat)
    ## ratio of preference
    rat <- spmat * (clsum[1] - clsum) /
        (clsum * sweep(-spmat, 2, spmat[1,], "+") + 1e-7)
    rat <- t(rat)
    ## finally construct species data
    spmat <- matrix(0, n, 3 * max(iclass))
    i <- seq(1, 3 * max(iclass), by = 3)
    spmat[,i] <- rat >= 0.8
    spmat[,i+1] <- rat >= 2
    spmat[,i+2] <- rat >= 6
    nc <- ncol(spmat)
    jdat <- numeric(ndat/2)
    k <- 1
    for (j in seq_len(n)) {
        for (i in seq_len(nc)) {
            if (spmat[j, i] > 0) {
                jdat[k] <- i
                k <- k + 1
            }
        }
        jdat[k] <- -1
        k <- k + 1
    }
    ## clsum up to every cluster and to same class level
    icl <- rep(seq_along(clsum), each = 3)
    clsum <- rep(clsum, each = 3)
    lev <- 2^levmax
    while(any(up <- icl < lev)) {
        clsum[up] <- sqrt(2) * clsum[up]
        icl[up] <- 2 * icl[up]
    }
    zidat <- Z$idat ## full data
    ## second and third pseudoquadrat level get double weight
    clsum <- c(1, 2, 2) * clsum
    Z <- .Fortran("class", nspec=as.integer(n),
                  icmax=as.integer(nc),
                  ndat=as.integer(length(jdat)),
                  0L, mmz=as.integer(MMZ), mms=as.integer(MMS), ix=Z$ix,
                  jnam=Z$jnam, iirow=Z$iirow, iaddr=which(jdat == -1),
                  indpot=seq_len(nc),
                  indord=Z$indord, izone=Z$izone, iy=Z$iy, jjcol=Z$jjcol,
                  jdat=as.integer(jdat), indsig=Z$indsig, ipict=Z$ipict,
                  x=Z$x, xx=Z$xx, rtot=Z$rtot, rrwt=as.double(spsum),
                  rowwgt=Z$rowwgt,
                  y=Z$y, yy=Z$yy, ctot=Z$ctot, ccwt=as.double(clsum),
                  colwgt=Z$colwgt,
                  jname1=Z$jname1, inflag=as.integer(inflag),
                  x3=Z$x3, x4=Z$x4, x5=Z$x5, lind=Z$lind,
                  inlevmax = Z$inlevmax, inmmin=Z$inmmin,
                  eig=double(2^levmax-1), indics = integer(inddim),
                  limpos = integer(2^levmax-1), isec=2L)
    species <- list(iclass = Z$jnam[seq_len(n)], eig = Z$eig,
                    labels = cnames)
    ## ordered index for quadrats and species: classorder() below
    quadrat$index <- classorder(quadrat$iclass, levmax)
    sindex <- classorder(species$iclass, levmax)
    irev <- 0
    rev <- .Fortran("revspec", nspec=as.integer(n), mm=as.integer(mm),
                    ndat=Z$ndat, ix=as.integer(quadrat$index),
                    iy=as.integer(sindex),
                    x=Z$x, y=Z$y, idat=zidat,
                    irev=as.integer(irev))$irev
    if (rev)
        sindex <- rev(sindex)
    species$index <- sindex
    ## out
    out <- list(call = match.call(), cutlevels = cutlevels,
                levelmax = levmax,
                nspecies = n, nquadrat = mm, idat = idat,
                quadrat = quadrat, species = species)
    class(out) <- "twinspan"
    out
}

### Internal function to sort items by their class number. For
### sorting, the class numbers must be elevated to the same level.

`classorder` <-
    function(class, levmax)
{
    ## all classes to the uppermost level numbers
    up <- 2^levmax
    while(any(small <- class < up)) {
        class[small] <- class[small] * 2
    }
    order(class)
}
