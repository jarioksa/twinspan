### sum of indicator scores for each quadrat

#' Indicator Scores for Quadrats
#'
#' Indicator scores are the sum of negative and positive pseudospecies
#' indicators on each quadrat.
#'
#' Function calculates the indicator scores that are used to split
#' quadrats at each division. These scores are the sum of positive and
#' negative indicator pseudospecies at each division. The indicator
#' pseudospecies can be seen in \code{\link{summary.twinspan}} and the
#' actual indicator values for each pseudospecies with
#' \code{\link{indvalues}}.
#'
#' Function has a \code{plot} method that display indicator scores
#' against first ordination axis using different colours for the
#' negative (left) and positive (right) group of the subdivision. The
#' \code{plot} adds a horizontal line that separates negative (left)
#' and positive (right) groups, and vertical lines for tentative
#' \dQuote{indifferent} zone for ordination split. If a quadrat is
#' deep in the negative or positive side of the ordination axis, it
#' will be classified by ordination against its indicator
#' scores. These quadrats appear as points with \dQuote{wrong} colour
#' for their indicator score, and they are said to be
#' \dQuote{misclassified} by indicators (these quadrats are reported
#' by function \code{\link{misclassified}}). \code{twinspan} uses
#' \dQuote{polished} ordination for final classification instead of
#' raw ordination axis. This is not accessible with the current code,
#' and therefore the vertical lines are tentative.
#'
#' @seealso \code{\link{summary.twinspan}} for divisions and indicator
#'     pseudospecies, \code{\link{indvalues}} for actual indicator
#'     values of pseudospecies, and \code{\link{misclassified}} for
#'     misclassified quadrats.
#'
#' @return Function returns an object of class \code{"indscores"}
#'     which is a list of following items:
#' \describe{
#'     \item{score}{Indicator score for each quadrat.}
#'     \item{classification}{Classes within this division.}
#'     \item{ordination}{First axis of correspondence analysis within
#'     this division.}
#'     \item{positivelimit}{Lower limit of indicator
#'     score for the positive (right) group.} }
#'
#' @examples
#' data(ahti)
#' tw <- twinspan(ahti)
#' summary(tw, level = 3, maxitems = 6)
#' indscores(tw, division = 8)
#' summary(indvalues(tw, division = 8))
#' misclassified(tw)
#' plot(indscores(tw, 1))
#' plot(indscores(tw, 3))
#'
#' @param object \code{twinspan} result object.
#' @param division Number code of \code{twinspan} division for quadrats.
#'
#' @importFrom stats cor
#' @importFrom vegan decorana scores
#'
#' @export
`indscores` <-
    function(object, division)
{
    if(object$quadrat$eig[division] == 0 ||
       length(object$quadrat$eig) < division)
        stop(division, " is not a division")
    ## indicator scores
    x <- twin2stack(object, subset = twingroup(object, division),
                    downweight = FALSE)
    inds <- object$quadrat$indicators[, division]
    inds <- inds[inds != 0]
    labs <- object$quadrat$indlabels[abs(inds)]
    x <- x[, labs, drop=FALSE]
    x <- sweep(x, 2, sign(inds), "*")
    sco <- rowSums(x)
    ## split of division
    level <- floor(log2(division))
    cl <- cut(object, level + 1)[twingroup(object, division)]
    ## first ordination axis
    ord <- decorana(twin2stack(object, subset = twingroup(object, division),
                               downweight = TRUE),
                    ira = 1, iresc = 0)
    ord <- drop(scores(ord, display = "sites", choices = 1))
    if (cor(ord, sco) < 0)
        ord <- -ord
    ## out
    structure(list(score = sco,
                   classification = cl,
                   ordination = ord,
                   positivelimit = object$quadrat$positivelimit[division]),
              class = c("indscores", "list"))
}

#' @rdname indscores
#'
#' @param x \code{indscores} result object.
#' @param col Colours for negative (left) and positive (right) groups.
#' @param \dots Other arguments passed to \code{plot}.
#'
#' @importFrom graphics abline
#'
#' @export
`plot.indscores` <-
    function(x, col = c("blue", "red"), ...)
{
    col <- col[x$classification - min(x$classification) + 1]
    cl <- sort(unique(x$classification))
    maxleft <- max(x$ordination[x$classification == cl[1]])
    minright <- min(x$ordination[x$classification == cl[2]])
    plot(score ~ ordination, data = x, col = col, ...)
    if (minright > maxleft)
        abline(v = (maxleft + minright) / 2)
    else {
        segments(maxleft, min(x$score) - 0.5, maxleft, x$positivelimit - 0.5)
        segments(maxleft, x$positivelimit - 0.5, maxleft, max(x$score) + 0.5,
                 lty = 2)
        segments(minright, min(x$score) - 0.5, minright, x$positivelimit - 0.5,
                 lty = 2)
        segments(minright, x$positivelimit - 0.5, minright, max(x$score + 0.5))
    }
    abline(h = x$positivelimit - 0.5)
}
