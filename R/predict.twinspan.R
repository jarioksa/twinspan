### predict classification using official indicator pseudospecies

#' Predict Class Membership of Quadrats
#'
#' Function predicts the class membership for each quadrat using the
#' reported indicator pseudospecies and limit for the indicator score
#' for the \dQuote{positive} (right) group.
#'
#' The \code{twinspan} classification is based on splicing polarized
#' ordination axis, and the reported indicator pseudospecies only
#' indicate the decisions in each division, and do not necessarily
#' give the same classification: The original classification cannot be
#' necessarily found when giving the original data as
#' \code{newdata}. In the original TWINSPAN this is called
#' misclassification.
#'
#' @seealso \code{\link{cut.twinspan}} gives the original
#'     classification, and \code{\link{misclassified}} analyses the
#'     differences of this and \code{predict}.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' predict(tw)
#' predict(tw, level=3)
#' ## misclassifications: predict and twinspan differ
#' sum(predict(tw) != cut(tw))
#' ## build model for 4/5 of data and predict for the removed 1/5
#' i <- rep(1:5, length = nrow(ahti))
#' i <- sample(i) # shuffle in random order
#' tw <- twinspan(ahti[i != 1,]) # remove i==1
#' predict(tw, newdata = ahti[i==1,])
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
        level <- 15
    inds <- object$quadrat$indicators
    poslim <- object$quadrat$positivelimit
    indlab <- object$quadrat$indlabels
    cuts <- object$cutlevels
    ## handle newdata
    if (missing(newdata))
        newdata <- twin2stack(object)
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
