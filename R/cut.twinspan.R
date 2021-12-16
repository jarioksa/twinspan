#' Return twinspan Classification at Given Level
#'
#' Function returns a vector of \code{twinspan} classes at given level
#' of hierarcy for quadrats or species.
#'
#' \code{\link{twinspan}} returns only the the classification at the
#' final level, but any upper level classes can be found by integer
#' divisions by 2. \code{\link{twinspan}} bases classification
#' principally on splitting polished ordination axis. Sometimes this
#' allocation is in conflict with indicator pseudospecies. This is
#' called \sQuote{misclassification} (see
#' \code{\link{misclassified}}). Function
#' \code{\link{predict.twinspan}} returns the similar classification
#' based on indicator pseudospecies (also for new data), and
#' \code{\link{misclassified}} analyses the differences of these
#' classifications.
#'
#' @seealso \code{\link{predict.twinspan}} gives similar classes, but
#'     based on indicator pseudospecies. \code{\link{cutree}} provides
#'     a similar functionality for \code{\link{hclust}} trees.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' cut(tw)
#' cut(tw, level=3)
#' cut(tw, what = "species")
#'
#' @param x \code{twinspan} result.
#' @param level Level of hierarchy for classification. If missing, the
#'     final level used in the object will be returned.
#' @param what Return either a \code{"quadrat"} or \code{"species"}
#'     classification vector.
#' @param \dots Other parameters (ignored).
#'
#' @return A vector of class numbers for the given level of hierarchy
#'     using Twinspan identifiers. For identifiers and levels, see
#'     \code{\link{summary.twinspan}},
#'     \code{\link{as.hclust.twinspan}}.
#'
#' @export
`cut.twinspan` <-
    function(x, level, what = c("quadrat", "species"), ...)
{
    what <- match.arg(what)
    if (missing(level))
        level <- x$levelmax
    clmax <- 2^(level+1) - 1
    cl <- x[[what]]$iclass
    while(any(big <- cl > clmax)) {
        ## mother class by integer division
        cl[big] <- cl[big] %/% 2
    }
    cl
}

## cut by group homogeneity as defined by within-group
## Chi-square. This could give similar clustering as Rolecek's
## modified twinspan which only splits the most heterogeneous group at
## each step instead of dichotomizing. However, we use the internal
## twinspan criterion of Chi-square instead of dissimilarities. There
## is no guarantee that the tree will be ordered, but let's hope so
## and as.dendrogram(x, height="chi") will show this.

#' @export
`chicut` <-
    function(x, ngroups)
{
    if (missing(ngroups))
        ngroups <- 1 # return one group if nothing is asked
    chi <- twintotalchi(x)
    ## latter half of chi have items that cannot be split
    k <- length(chi) %/% 2L + 1L
    chi[k:length(chi)] <- 0
    ## terminal nodes (leaves) cannot be split
    chi[x$quadrat$eig <= 0] <- 0
    ix <- order(chi, decreasing = TRUE) # order by heterogeneity
    clim <- 2^(0:x$levelmax) - 1L
    class <- rep(1, x$nquadrat)
    for(i in seq_len(ngroups - 1)) {
        lev <- max(which(ix[i] > clim))
        id <- 2L * ix[i]
        class[cut(x, level=lev) == id] <- id
        class[cut(x, level=lev) == id+1L] <- id + 1L
    }
    class
}
