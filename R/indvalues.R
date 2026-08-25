### evaluate indicator scores for pseudospecies in a division
### similarly as in twinspan subroutine indsco

#' Indicator Scores of Pseudospecies in a Division
#'
#' Function returns the indicator scores that \code{\link{twinspan}}
#' uses to select indicator pseudospecies.
#'
#' Indicator scores are differences of relative frequencies of
#' pseudospecies in the posivite (right) and negative (left)
#' branches. Perfect indicator pseudopecies has value 1 or -1 and
#' occurs on every quadrat of one branch and in none of the other
#' branch. Indicator pseudospecies of division are selected by the
#' indicator scores, but not all obviously best species occur as
#' indicators. One basic species can be listed only once using its
#' best pseudospecies, and redundant pseudospecies are not listed.
#'
#' @seealso Function \code{\link{summary.twinspan}} lists the used
#'     indicator pseudospecies in each division.
#'
#' @examples
#' data(ahti)
#' tw <- twinspan(ahti)
#' summary(indvalues(tw, 1))
#' summary(indvalues(tw, 2))
#' summary(indvalues(tw, 3))
#'
#' @param object \code{twinspan} result object.
#' @param division Number code of \code{twinspan} division for quadrats.
#'
#' @export
`indvalues` <-
    function(object, division)
{
    if(object$quadrat$eig[division] == 0 ||
       length(object$quadrat$eig) < division)
        stop(division, " is not a division")
    x <- twin2stack(object, subset = twingroup(object, division),
                    downweight = FALSE)
    level <- floor(log2(division))
    cl <- cut(object, level + 1)[twingroup(object, division)]
    freq <- by(x, cl, colMeans)
    inds <- freq[[2]] - freq[[1]]
    class(inds) <- "indvalues"
    inds
}

#' @rdname indvalues
#' @param value Minimum absolute value of returned indicator scores.
#' @param \dots Other arguments passed to the function (ignored).
#' @export
`summary.indvalues` <-
    function(object, value = 0.25, ...)
{
    cat("\nBest Indicator Pseudospecies:\n")
    object <- object[abs(object) >= value]
    object[rev(order(abs(object)))]
}
