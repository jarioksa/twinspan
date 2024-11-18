#' Transform Data for Correspondence Analysis like twinspan
#'
#' Function transforms data so that Correspondence Analysis gives the
#' same result as in \code{\link{twinspan}} divisions.
#'
#' In \code{\link{twinspan}}, quantitative species data are split into
#' binary (0/1) pseudospecies by \code{cutlevels}. All these
#' pseudospecies are stacked as columns in a new data set. Rare
#' pseudospecies that occur at lower frequency than 0.2 are
#' \code{\link[vegan]{downweight}}ed within \code{twinspan}. This
#' reduces the weight of rare species or rare abundance levels in
#' correspondence analysis, but downweighting is optional in this
#' function.
#'
#' When the downweighted data are analysed with correspondence
#' analysis (e.g., \code{\link[vegan]{cca}},
#' \code{\link[vegan]{decorana}} with option \code{ira=1}), these will
#' give the same first eigenvalue and ordination as in
#' \code{\link{twinspan}}. When a \code{subset} of a
#' \code{\link{twinspan}} class is used, correspondence analysis of
#' subdivision of the class can be obtained.
#'
#' @seealso \code{\link[vegan]{downweight}} in \CRANpkg{vegan}: this
#'     function is often used with Detrended Correspondence Analysis
#'     (\code{\link[vegan]{decorana}}). However, the implementation is
#'     slightly different in TWINSPAN, and weights differ
#'     slightly. Function \code{\link{twin2stack}} extracts similar
#'     data from a \code{\link{twinspan}} result object.
#'
#' @examples
#'
#' data(ahti)
#' tahti <- twinsform(ahti)
#' colnames(tahti)
#' ## needs vegan for correspondence analysis
#' if (suppressPackageStartupMessages(require("vegan"))) {
#' decorana(tahti, ira=1)
#' }
#' ## similar first eigenvalue
#' eigenvals(twinspan(ahti))
#'
#' @return A stacked matrix of optionally downweighted pseudospecies.
#'
#' @param x Input (community) data.
#' @param cutlevels Cut levels used to split quantitative data into
#'     binary pseudospecies.
#' @param subset Logical vector or indices that select a subset of
#'     quadrats (sampling units).
#' @param downweight Downweight result similarly as in
#'     \code{\link[vegan]{decorana}}. Downweighting is needed to
#'     replicate the process in \code{twinspan}, but it can be left
#'     out when we only want to have a stacked data set for other
#'     uses.
#'
#' @export
`twinsform` <-
    function(x, cutlevels = c(0,2,5,10,20), subset, downweight = TRUE)
{
    x <- as.matrix(x)
    ## take a subset
    if (!missing(subset))
        x <- x[subset,, drop = FALSE]
    ## Twinspan stacks binary matrices of pseudospecies and
    ## downweights them
    nlev <- length(cutlevels)
    nm <- colnames(x)
    if (cutlevels[1] <= 0)
        ax <- x > cutlevels[1]
    else
        ax <- x >= cutlevels[1]
    for (k in 2L:nlev)
        ax <- cbind(ax, x >= cutlevels[k])
    colnames(ax) <- paste0(nm, rep(seq_len(nlev), each=length(nm)))
    ax <- ax[, colSums(ax) > 0, drop = FALSE]
    ## vegan::downweight implements decorana downweighting, but
    ## TWINSPAN uses slightly different (and more approximate form).
    if (downweight) {
        ax <- downweight(ax)
    }
    ax
}

### internal support function for downweighting the TWINSPAN way
### (which is different than in DECORANA, and pretty crude way of
### doing things).

## @param x input data matrix
## @param rw prior row weights

`downweight` <-
    function(x, rw = rep(1, NROW(x)))
{
    cs <- colSums(rw * x)
    v <- rep(1, ncol(x))
    lim <- sum(rw)/5
    down <- cs < lim
    v[down] <- (cs/lim)[down] * 0.99 + 0.01
    x <- sweep(x, 2, v, "*")
    attr(x, "v") <- v
    attr(x, "fraction") <- 5 # for vegan::downweight compatibility
    x
}
### A non-exported internal function to construct stacked community
### matrix from internal data returned by twinspan. The output is
### similar as with twinsform(x, downweight = FALSE).

#' @rdname twin2mat
#'
#' @title Extract Transformed Input Data from twinspan Result
#'
#' @description
#'
#' Functions extract the data \code{\link{twinspan}} used in its
#' analysis, and allow reproducing the internal ordination and
#' inspecting the \code{twinspan} divisions.
#'
#' @details
#'
#' Function \code{twin2stack} extracts the pseudospecies matrix, where
#' columns are pseudospecies with their cutlevels. This is similar to
#' the file generated with \code{\link{twinsform}}. The default is to
#' return a binary matrix, where data entries are eiter \eqn{0} or
#' \eqn{1}. Alternatively, it is possible to extract a subset of data
#' with downweighting allowing scrutiny of \code{\link{twinspan}}
#' divisions. When downweighted data are ordinated with correspondence
#' analysis (such as \CRANpkg{vegan} functions
#' \code{\link[vegan]{cca}}, \code{\link[vegan]{decorana}} with
#' \code{ira=1}) the first eigenvalue will match the eigenvalue in
#' \code{\link{twinspan}}, and when a division is used as a
#' \code{subset}, its eigenvalue will match with \code{twinspan}.
#'
#' Function \code{twin2mat} extracts data file with pseudospecies
#' transformation. Columns are original species, and entries are
#' abundances after pseudospecies transformation. This is similar as
#' the output from \CRANpkg{vegan} function
#' \code{\link[vegan]{coverscale}} with similar cut levels and
#' argument \code{character=FALSE}. These data were not analysed in
#' \code{\link{twinspan}}, but these are the data tabulated with
#' \code{\link{twintable}}.
#'
#' Function \code{twin2specstack} returns similar data as used in
#' species classification in \code{\link{twinspan}}. In this matrix,
#' species are rows and columns are \dQuote{pseudocluster}
#' preferences. The preference of each species in each terminal group
#' and internal division is estimated as proportion of its abundance
#' (in pseudospecies scale, see \code{twin2mat}) in the group and all
#' data. If this proportion is 0.8 or higher, species is regarded as
#' present at pseudocluster value 1, if it is 2 or higher at value 2,
#' and if it is 6 or higher at value 3. With \code{downweight=FALSE}
#' these data are returned. The columns are named by their division or
#' cluster number followed \code{a}, \code{b} and \code{c} for
#' pseudocluster levels (and including zero columns). In default, the
#' pseudocluster values are still downweighted using species
#' frequencies as weights, and then rows are weighted by species
#' frequencies and columns by their totals extended to the same lowest
#' level of classification, giving two times higher weight to higher
#' \dQuote{pseudocluster} levels \code{b} and \code{c}.  When
#' ordinated with correspondence analysis (\code{\link[vegan]{cca}},
#' \code{\link[vegan]{decorana}} with \code{ira=1}) this gives similar
#' eigenvalue for the first axis as in \code{\link{twinspan}}, and
#' when a division is used as a \code{subset}, similar eigenvalue as
#' in that division.
#'
#' @seealso For original data set instead of \code{\link{twinspan}}
#'     result, functions \code{\link{twinsform}} and
#'     \code{\link[vegan]{coverscale}} are analogous to
#'     \code{twin2stack} and \code{twin2mat}.
#'
#' @examples
#'
#' data(ahti)
#' dim(ahti)
#' range(ahti)
#' tw <- twinspan(ahti)
#' x <- twin2mat(tw)
#' dim(x)
#' range(x)
#' colnames(x)
#' x <- twin2stack(tw)
#' dim(x)
#' range(x)
#' colnames(x)
#' ## Inspect group 4
#' x <- twin2stack(tw, subset = twingroup(tw, 4), downweight = TRUE)
#' ## need vegan for correspondence analysis
#' if (suppressPackageStartupMessages(require("vegan"))) {
#' cca(x)
#' }
#' ## species classification
#' x <- twin2specstack(tw)
#' if (suppressPackageStartupMessages(require("vegan"))) {
#' cca(x)
#' }
#'
#' @param x \code{\link{twinspan}} result object.
#' @param subset Select a subset of quadrats (\code{twin2stack}) or
#'     species (\code{twin2specstack}).
#' @param downweight Downweight infrequent pseudospecies.
#'
#' @export
`twin2stack` <-
    function(x, subset, downweight = FALSE)
{
    nc <- length(x$quadrat$indlabels)
    nr <- x$nquadrat
    idat <- x$idat
    out <- matrix(0L, nr, nc)
    dimnames(out) <- list(x$quadrat$labels, x$quadrat$indlabels)
    i <- 1
    ## Internally twinspan stores data as a vector of pseudospecies
    ## indices (each numerically with abundance 1) separating quadrats
    ## (SUs) by -1.
    for(j in seq_along(idat)) {
        ## new quadrat?
        if(idat[j] == -1) {
            i <- i + 1
            next
        }
        ## pseudospecies is present with abundance 1
        out[i, idat[j]] <- 1L
    }
    if (!missing(subset)) {
        out <- out[subset,, drop = FALSE]
        cs <- colSums(out)
        if (any(cs==0))
            out <- out[, cs > 0, drop = FALSE]
    }
    if (downweight)
        out <- downweight(out)
    out
}

### A non-exported internal function to construct pseudospecies
### community matrix from internal data returned by twinspan. The
### numeric values are pseudo-species transformed values. The output
### is similar as with vegan::coverscale(x, "Hill",
### character=FALSE). The function uses the same internal data as
### twin2stack(), but adds binary pseudospecies to transformed values.

#' @rdname twin2mat
#' @export
`twin2mat` <-
    function(x)
{
    nc <- x$nspecies
    nr <- x$nquadrat
    idat <- x$idat
    out <- matrix(0L, nr, nc)
    dimnames(out) <- list(x$quadrat$labels, x$species$labels)
    ## Indices of pseudospecies > 1 are stored after the first level
    ## of pseudospecies, for which the index matches the original
    ## data. We update species index only when we meet such an
    ## original index, and add values to matrix at every step.
    i <- 1
    for(k in seq_along(idat)) {
        ## new quadrat?
        if (idat[k] == -1) {
            i <- i + 1
            next
        }
        ## new species?
        if (idat[k] <= nc) {
            j <- idat[k]
        }
        out[i,j] <- out[i,j] + 1L
    }
    out
}

### Reconstruct data for getting species clustering in TWINSPAN. In
### these data, species are the observations (rows), and clusters &
### divisions (internal nodes) are "pseudoquadrats". For each division
### & cluster, we evaluate the preference of species for that
### cluster. If it occurs at >= 0.8 frequency to the average (little
### less than usually), it is regarded as present at
### "pseudoquadrat1". If it is 2x times more frequent, it is also
### "pseudoquadrat2", and 6x preferentials are "pseudoquadrat3". So
### each cluster/division is divided into three pseudoquadrats for
### each species, solely spaced on species frequencies. Then these
### data are downweighted with row weights of species frequency, and
### for CA both the species (rows) and the pseudoquadrats are weighted
### by the cluster total occurrences of species, and the
### "pseudoquadrat2" and "pseudoquadrat3" get double weight to
### "pseudoquadrat1".

#' @rdname twin2mat
#' @export
`twin2specstack` <-
    function(x, subset, downweight = TRUE)
{
    cl <- x$quadrat$iclass
    nsp <- x$nspecies
    mat <- twin2mat(x)  # pseudospecies data
    ## collect counts for all groups and their mother divisions
    clmax <- max(cl)
    tot <- numeric(clmax)
    ## terminal groups
    rs <- rowSums(mat)
    for (k in seq_along(cl)) {
        tot[cl[k]] <- tot[cl[k]] + rs[k]
    }
    ## mother divisions
    for (k in clmax:2) {
        tot[k %/% 2] <- tot[k %/% 2] + tot[k]
    }
    ## same for all species: first terminal groups
    totj <- matrix(0, clmax, nsp)
    dimnames(totj) <- list(seq_len(clmax),x$species$labels)
    for (k in seq_along(cl)) {
        totj[cl[k],] <- totj[cl[k],] + mat[k,]
    }
    ## mother divisions
    for (k in clmax:2) {
        totj[k %/% 2,] <- totj[k %/% 2,] + totj[k,]
    }
    ## ratio of preference for a class
    rat <- totj * (tot[1] - tot) /
        (tot * sweep(-totj, 2, totj[1,], "+") + 1e-7)
    ## we transpose: groups as rows nicer for vector * matrix
    ## calculations, but now we get species as rows
    rat <- t(rat)
    ## binary species matrix: TWINSPAN uses arbitrary constants 0.8, 2
    ## and 6 for pseudovalues
    smat <- cbind(rat >= 0.8, rat >= 2, rat >= 6) + 0
    ## rename pseudoquadrat levels as a, b, c
    colnames(smat) <- paste0(colnames(smat), rep(letters[1:3], each=clmax))
    ## get out if downweight = FALSE
    if (!downweight)
        return(smat)
    ## row (species) species weights are species frequencies and
    ## column (cluster) weights are cluster totals elevated to the
    ## same level of division
    rwt <- colSums(mat > 0)
    cwt <- tot
    icl <- seq_len(clmax)
    lev <- 2^(x$levelmax)
    while(any(up <- icl < lev)) {
       cwt[up] <- sqrt(2) * cwt[up] # TWINSPAN multiplier is 1.414
       icl[up] <- 2 * icl[up]
    }
    ## First (a) pseudoquadrat level gets this weight, upper levels
    ## (b, c) twice the weight
    cwt <- c(cwt, 2*cwt, 2*cwt)
    ## subset of species
    if (!missing(subset)) {
        smat <- smat[subset,, drop=FALSE]
        rwt <- rwt[subset]
    }
    ## weighted downweighting
    smat <- downweight(smat, rw = rwt)
    ## row and column weighting for CA
    smat <- rwt * sweep(smat, 2, cwt, "*")
    smat <- smat[, colSums(smat) > 0, drop = FALSE]
    smat
}
