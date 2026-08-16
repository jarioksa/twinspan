### twinspan interface to vegan::tabasco

#' Display twinspan Classification as a Heatmap
#'
#' Display \code{twinspan} classification as a \code{link{heatmap}}
#' with quadrat and species dendrograms using \pkg{vegan} function
#' \code{\link[vegan]{tabasco}}.
#'
#' @seealso \code{\link{twintable}}, \code{\link{image.twinspan}},
#'     \code{\link[vegan]{tabasco}}.

#' @inheritParams twintable
#' @param height Use either division levels (\code{"level"}), total
#'     Chi-squares (\code{"chi"}) or eigenvalues of first axis of
#'     division (\code{"eigen"}) as dendrogram heights.
#' @param Rowv,Colv Reorder dendrograms. If \code{FALSE} use original
#'     \code{twinspan} order. If \code{TRUE} order by the first axis
#'     of CCA. Alternatively, the arguments can be vectors that are
#'     used to reorder the dendrogram.
#' @param \dots Other arguments passed to \code{\link[vegan]{tabasco}}.
#'
#' @importFrom vegan tabasco
#'
#' @export
`twinheat` <-
    function(object, maxspp, goodspecies, subset, height = "level",
             Rowv = FALSE, Colv = FALSE, ...)
{
    if (!inherits(object, "twinspan"))
        stop("function can be used only with 'twinspan' result object")
    x <- twin2mat(object)
    if (!missing(subset)) {
        if (!is.logical(subset))
            subset <- sort(subset) # we do not want to reorder by subset
        x <- x[subset,, drop=FALSE]
    }
    present <- colSums(x) > 0
    skept <- logical(length(present))
    if (!missing(goodspecies)) {
        if (!missing(maxspp))
            stop("maxspp and goodspecies cannot be defined together")
        k <- twinspan:::goodspec(object, goodspecies)
        skept[k] <- TRUE
        skept <- skept & present
    }
    else if (!missing(maxspp) && sum(present) > maxspp) {
        sptot <- colSums(x)
        maxlim <- sort(sptot, decreasing = TRUE)[maxspp]
        skept <- sptot >= maxlim & present
    } else {
        skept <- present
    }
    sptree <- as.dendrogram(object, what = "species", subset = skept,
                            height = height)
    qtree <- as.dendrogram(object, subset = subset, height = height)
    tabasco(x, use = qtree, sp.ind = sptree, Rowv = Rowv, Colv = Colv, ...)
}
