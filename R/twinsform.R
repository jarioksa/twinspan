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
        x <- x[subset, ]
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

`downweight` <-
    function(x)
{
    cs <- colSums(x)
    v <- rep(1, ncol(x))
    lim <- nrow(x)/5
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
#' divisions.
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
#' ## Inspect division 4
#' x <- twin2stack(tw, subset = cut(tw, 2) == 4, downweight = TRUE)
#' ## need vegan for correspondence analysis
#' if (suppressPackageStartupMessages(require("vegan"))) {
#' cca(x)
#' }
#'
#' @param x \code{\link{twinspan}} result object.
#' @param subset Select a subset of quadrats.
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
