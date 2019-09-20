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

`twin2stack` <-
    function(x)
{
    nc <- length(x$quadrat$indlabels)
    nr <- x$nquadrat
    idat <- x$idat
    out <- matrix(0, nr, nc)
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
