#' Transform Data for Correspondence Analysis like twinspan
#'
#' Function transforms data so that Correspondence Analysis gives the
#' same result as in \code{\link{twinspan}}.
#'
#' In \code{\link{twinspan}}, quantitative species data are split into
#' binary (0/1) pseudospecies by \code{cutlevels}. All these
#' pseudospecies are stacked as columns in a new data set. Rare
#' pseudospecies that occur at lower frequency than 0.2 are
#' \code{\link[vegan]{downweight}}ed to reduce the importance of rare
#' species or rare abundance levels in correspondence analysis. When
#' analysed with correspondence analysis (e.g.,
#' \code{\link[vegan]{cca}}, \code{\link[vegan]{decorana}} with option
#' \code{ira=1}), this will give the same eigenvalue and ordination as
#' used in \code{\link{twinspan}}. When a \code{subset} of a
#' \code{\link{twinspan}} class is used, correspondence analysis of
#' subdivision of the class can be obtained. The results are not
#' absolutely similar, mainly because correspondence analysis is more
#' approximate in \code{\link{twinspan}} than in dedicated
#' functions. Moreover, \code{\link{twinspan}} adds a small value
#' (0.01) to all abundance classes, and does not use weights smaller
#' than 0.01 which can influence results when the number of quadrats
#' is above 100.
#'
#' @examples
#'
#' data(ahti)
#' tahti <- twinsform(ahti)
#' colnames(tahti)
#' ## needs vegan for correspondence analysis
#' if (require(vegan)) {
#' ## similar first eigenvalue
#' decorana(tahti, ira=1)
#' eigenvals(twinspan(ahti))
#' }
#'
#' @return A stacked matrix of downweighted pseudospecies.
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
#' @importFrom vegan downweight
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
    ax <- ax[, colSums(ax) > 0]
    if (downweight)
        ax <- downweight(ax)
    ax
}

### A non-exported internal function to construct stacked community
### matrix from internal data returned by twinspan. The output is
### similar as with twinsform(x, downweight = FALSE).

#' @rdname twin2mat
#'
#' @title Extract Transformed Input Data from TWINSPAN result
#'
#' @description
#'
#' Function \code{twin2stack} extracts the binary pseudospecies file,
#' where columns are pseudospecies with their cutlevels, and entries
#' are either \code{0} or \code{1}. This is similar to the file
#' generated with \code{\link{twinsform}} with argument
#' \code{downweight=FALSE}. Function \code{twin2mat} extracts data
#' file with pseudospecies transformation. Columns are original
#' species, and entries are abundances after pseudospecies
#' transformation. This is similar as the output from \CRANpkg{vegan}
#' function \code{\link[vegan]{coverscale}} with similar cut levels
#' and argument \code{character=FALSE}.
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
#'
#' @param x \code{\link{twinspan}} result object.
#'
#' @export
`twin2stack` <-
    function(x)
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
