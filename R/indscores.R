### sum of indicator scores for each quadrat

#' Indicator Scores for Quadrats
#'
#' Indicator scores are the sum of negative and positive pseudospecies
#' indicators on each quadrat.
#'
#' @inheritParams indvalues
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
              clas = "indscores")
}

#' @rdname indscores
#'
#' @param x \code{indscores} result object.
#' @param col Colours for negative (left) and positive (right) groups.
#'
#' @importFrom graphics abline
#'
#' @export
`plot.indscores` <-
    function(x, col = c("red", "blue"), ...)
{
    col <- col[x$classification - min(x$classification) + 1]
    plot(score ~ ordination, data = x, col = col, ...)
    abline(v = min(x$ordination)/5)
    abline(v = max(x$ordination)/5)
    abline(h = x$positivelimit - 0.5)
}
