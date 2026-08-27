### predict classification using official indicator pseudospecies

#' Predict Class Membership of Quadrats
#'
#' Function predicts the class membership for each quadrat using the
#' reported indicator pseudospecies and limit for the indicator score
#' for the \dQuote{positive} (right) group.
#'
#' The \code{twinspan} classification is based on splitting ordination
#' axis, and the indicator pseudospecies are used for decision in the
#' \dQuote{indifferent zone} on the border of two classes. Points
#' further towards the ends are classified by their position on the
#' ordination axis. The \code{predict} function uses indicator scores
#' and the predicted classification differs for the points that were
#' primarily classified by their position on the ordination
#' axis. Function \code{\link{misclassified}} reports these points,
#' and the classification can be visualized with the \code{plot} of
#' the \code{\link{indscores}}.
#'
#' Function can use \code{newdata} to predict classification for
#' quadrats that were not included in the \code{twinspan}
#' analysis. However, the method is not very good in extrapolating
#' classification with different set of species. It will happily
#' predict a class even with \code{newdata} that shares no species
#' with the \code{twinspan} result: it will follow the decision path
#' where every indicator value is 0. Function \code{dumpclass} will
#' report what would be the predicted class when \code{newdata} has no
#' indicator species. The members of that class should be carefully
#' inspected, and the result should be treated with caution if that
#' class is over-represented in the prediction.
#'
#' @seealso \code{\link{cut.twinspan}} gives the original
#'     classification, and \code{\link{misclassified}} analyses the
#'     differences of this and \code{predict}. Function
#'     \code{\link{indscores}} returns the indicator scores for
#'     quadrats in any division, and its \code{plot} function can be
#'     used for visualization of predictions.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' predict(tw)
#' dumpclass(tw)
#' predict(tw, level=3)
#' dumpclass(tw, level=3)
#' ## misclassifications: predict and twinspan differ
#' sum(predict(tw) != cut(tw))
#' ## build model for 4/5 of data and predict for the removed 1/5
#' i <- rep(1:5, length = nrow(ahti))
#' i <- sample(i) # shuffle in random order
#' tw <- twinspan(ahti[i != 1,]) # remove partition i==1
#' predict(tw, newdata = ahti[i==1,])
#' dumpclass(tw)
#'
#' @param object \code{twinspan} result object.
#' @param newdata Data used in prediction. The species will be matched
#'     by their names, and the pseudospecies are based on the
#'     \code{cutlevels} used in the original \code{twinspan} model.
#' @param level Level of hierarchy of classification. If missing, the
#'     prediction is made to the highest level of classification.
#' @param binname Use binary labels instead of decimal class numbers.
#' @param \dots Other parameters passed to the function (ignored).
#'
#' @importFrom stats predict
#'
#' @export
`predict.twinspan` <-
    function(object, newdata, level, binname = FALSE, ...)
{
    if (missing(level))
        level <- object$levelmax
    inds <- object$quadrat$indicators
    poslim <- object$quadrat$positivelimit
    indlab <- object$quadrat$indlabels
    cuts <- object$cutlevels
    ## handle newdata
    if (missing(newdata))
        newdata <- twin2stack(object, downweight = FALSE)
    else
        newdata <- twinsform(newdata, cuts, downweight = FALSE)
    ## one quadrat is turned into 1-column matrix: transpose
    if (NCOL(newdata) == 1)
        newdata <- t(newdata)
    m <- colnames(newdata) %in% indlab
    newdata <- newdata[,m, drop=FALSE]
    pred <- numeric(nrow(newdata))
    prow <- numeric(length(indlab))
    names(prow) <- indlab
    ## cycles over rows of newdata
    for(i in seq_len(nrow(newdata))) {
        prow[] <- 0
        prow[colnames(newdata)] <- newdata[i,]
        k <- 1
        for (lvl in seq_len(level)) {
            ind <- inds[,k]
            if (all(ind == 0)) break
            ind <- ind[ind != 0]
            score <- sum(sign(ind) * prow[abs(ind)])
            if (score < poslim[k]) {
                k <- 2*k
            } else {
                k <- 2*k + 1
            }
            pred[i] <- k
            if (k > ncol(inds)) break
        }
    }
    if (binname)
        pred <- sapply(pred, class2bin)
    pred
}

### predict will give a class even when there are no shared species
### between twinspan model and newdata: it will be the class when
### indicator score is 0 at every division. This function tells what
### would that class be.

#' @rdname predict.twinspan
#' @export
`dumpclass` <-
    function(object, level, binname = FALSE)
{
    if (missing(level))
        level <- object$levelmax
    ndiv <- 2^level - 1
    k <- 1
    repeat{
        if (k > ndiv || all(object$quadrat$indicators[,k] == 0))
            break
        if (0 < object$quadrat$positivelimit[k])
            k <- 2 * k
        else
            k <- 2 * k + 1
    }
    if (binname)
        k <- class2bin(k)
    k
}
