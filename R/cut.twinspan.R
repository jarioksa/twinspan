#' Return twinspan Classification at Given Level
#'
#' Returns a vector of \code{twinspan} classes at a given level of
#' hierarchy or classes respecting group heterogeneity for quadrats or
#' species.
#'
#' @encoding UTF-8
#'
#' @details
#'
#' \code{\link{twinspan}} returns only the the classification at the
#' final level, but any upper level classes can be found by integer
#' divisions by 2. Function \code{cut} returns a vector class id
#' numbers for a given level of classification. Utility function
#' \code{twingroup} returns a logical vector that is \code{TRUE} for
#' items belonging to a certain group at any level. It can be more
#' practical in subsetting data.
#'
#' Function \code{cuth} cuts the classification by class heterogeneity
#' instead of level, and can be used to implement the modified method
#' of Roleček et al. (2009). The groups are formed with decreasing
#' heterogeneity but respecting the hierarchy. Total chi-square (also
#' known as inertia) is used as the criterion of heterogeneity. The
#' criterion is calculated with \code{\link{twintotalchi}} and the
#' criterion is based on the same data matrix as internally used in
#' \code{\link{twinspan}}. The function can also be used for species
#' classification, also with the internally used modified species
#' matrix.
#'
#' @references
#'
#' Roleček, J, Tichý, L., Zelený, D. & Chytrý, M. (2009). Modified
#' TWINSPAN classification in which the hierarchy respects cluster
#' heterogeneity. \emph{J Veg Sci} 20: 596--602.
#'
#' @seealso \code{\link{predict.twinspan}} gives similar classes, but
#'     based on indicator pseudospecies. \code{\link{cutree}} provides
#'     a similar functionality for \code{\link{hclust}}
#'     trees. Function \code{\link{as.hclust.twinspan}} generates
#'     corresponding tree presentation, and
#'     \code{\link{plot.twinspan}} will print that tree labelling
#'     internal nodes (divisions).
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' cut(tw)
#' ## traditional twinspan classification by level of hierarchy
#' cut(tw, level=3)
#' cut(tw, what = "species")
#' ## number of groups as with level=3, but by group heterogeneity
#' cuth(tw, ngroups = 8)
#'
#' @param x \code{twinspan} result.
#' @param level Level of hierarchy for classification. If missing, the
#'     final level used in the object will be returned.
#' @param what Return either a \code{"quadrat"} or \code{"species"}
#'     classification vector.
#' @param binname Use binary label for classes instead of decimal number.
#' @param \dots Other parameters (ignored).
#'
#' @return A vector of class numbers for the given level of hierarchy
#'     using Twinspan identifiers. For identifiers and levels, see
#'     \code{\link{summary.twinspan}},
#'     \code{\link{as.hclust.twinspan}}.
#'
#' @export
`cut.twinspan` <-
    function(x, level, what = c("quadrat", "species"), binname = FALSE, ...)
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
    if (binname)
        cl <- sapply(cl, class2bin)
    cl
}

## utility functin to return logical vector that is TRUE for items in
## the given group.

#' @rdname cut.twinspan
#'
#' @param group Group id number.
#'
#' @export
`twingroup` <-
    function(x, group, what = c("quadrat", "species"))
{
    what <- match.arg(what)
    level <- floor(log2(group))
    cut(x, level, what) == group
}

## cut by group homogeneity as defined by within-group
## Chi-square. This could give similar clustering as Rolecek's
## modified twinspan which only splits the most heterogeneous group at
## each step instead of dichotomizing. There is no guarantee that the
## tree will be ordered, but the code should handle this.

#' @rdname cut.twinspan
#'
#' @param ngroups Number of groups.
#'
#' @export
`cuth` <-
    function(x, ngroups, what = c("quadrat", "species"), binname = FALSE)
{
    what <- match.arg(what)
    if (missing(ngroups))
        ngroups <- 1 # return one group if nothing is asked
    chi <- twintotalchi(x, what = what)
    ## latter half of chi have items that cannot be split
    k <- length(chi) %/% 2L + 1L
    chi[k:length(chi)] <- 0
    ## terminal nodes (leaves) cannot be split
    chi[x[[what]]$eig <= 0] <- 0
    ix <- order(chi, decreasing = TRUE) # order by heterogeneity
    ix <- ix[seq_len(sum(chi > 0))]
    ix <- fixTreeReversal(ix)
    clim <- 2^(0:x$levelmax) - 1L
    nobs <- switch(what,
                   "quadrat" = x$nquadrat,
                   "species" = x$nspecies)
    class <- rep(1, nobs)
    for(i in seq_len(ngroups - 1)) {
        lev <- max(which(ix[i] > clim))
        id <- 2L * ix[i]
        class[cut(x, what = what, level=lev) == id] <- id
        class[cut(x, what = what, level=lev) == id+1L] <- id + 1L
    }
    if (binname)
        class <- sapply(class, class2bin)
    class
}
