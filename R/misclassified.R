### report "misclassified" observations
#'
#' Report Misclassified Quadrats
#'
#' \code{\link{twinspan}} bases its quadrat classification primarily
#' on ordination axis. In some cases this is in conflict with the
#' classification derived from indicator pseudospecies. This function
#' identifies these cases.
#'
#' The function compares the final \code{\link{twinspan}}
#' classification (from \code{\link{cut.twinspan}}) and the
#' classification predicted from indicator pseudospecies (from
#' \code{\link{predict.twinspan}}). If these two differ, the quadrat
#' is \dQuote{misclassified}. With \code{verbose=FALSE}, the function
#' returns a logical vector that is \code{TRUE} for misclassified
#' quadrats. In default, it also returns the names, twinspan classes
#' and predicted classes of misclassified quadrats, and the division
#' where the misclassification occurred and classifications
#' diverged. The divisions and their numbers can be seen in
#' \code{\link{summary.twinspan}} and \code{\link{plot.twinspan}}.
#'
#' @return
#'
#' If \code{verbose=TRUE}, the function returns an object of class
#' \code{"misclassified"} with following elements.
#' \describe{
#' \item{index}{Index of misclassified quadrats.}
#' \item{labels}{Labels (names) of quadrats.}
#' \item{class}{Final classification from \code{\link{twinspan}}.}
#' \item{predicted}{Final classification from \code{\link{predict.twinspan}}.}
#' \item{division}{The division where the misclassification occurred.}
#' }
#'
#' With \code{verbose=FALSE}, the function returns a logical vector
#' that is \code{TRUE} for misclassified quadrats.
#'
#' @seealso The basic functions are \code{\link{cut.twinspan}} and
#'     \code{\link{predict.twinspan}}. You can see the division
#'     numbers with \code{\link{summary.twinspan}} (with indicator
#'     pseudospecies) and in \code{\link{plot.twinspan}}.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' misclassified(tw)
#' ## see the ID numbers of divisions
#' plot(tw)
#' ## only look at misclassifications at first two levels
#' misclassified(tw, level = 2)
#'
#' @param x \code{\link{twinspan}} result object.
#' @param level Only consider misclassification down to this level.
#' @param verbose If \code{FALSE}, returns only the index of
#'     misclassified quadrats.
#' @param binname Use binary labels instead of decimal numbers for
#'     classes and divergent division in verbose object.
#'
#' @export
`misclassified` <-
    function(x, level, verbose = TRUE, binname = FALSE)
{
    if (missing(level))
        level <- x$levelmax
    class <- cut(x, level = level)
    pred <- predict(x, level = level)
    # misclassified cases
    ind <- class != pred
    if (!verbose)
        return(ind)
    ind <- which(ind)
    ## take only misclassified cases
    labs <- x$quadrat$labels[ind]
    class <- class[ind]
    pred <- pred[ind]
    ## see the level of divergence
    level <- x$levelmax
    div <- numeric(length(ind))
    cl <- class
    pr <- pred
    for(lev in seq(level-1, 0)) {
        lim <- 2^(lev+1) - 1
        pr[pr > lim] <- pr[pr > lim] %/% 2L
        cl[cl > lim] <- cl[cl > lim] %/% 2L
        ## prediction and class match for the first time
        k <- pr == cl & div == 0
        if (any(k))
            div[k] <- pr[k]
        ## all done?
        if (all(div > 0))
            break
    }
    if (binname) {
        class <- sapply(class, class2bin)
        pred <- sapply(pred, class2bin)
        div <- sapply(div, class2bin)
    }
    out <- list(index = ind, labels = labs, class = class, predicted = pred,
                division = div)
    class(out) <- "misclassified"
    out
}

#' @export
`print.misclassified` <-
    function(x, ...)
{
    nmiss <- length(x$index)
    if (nmiss == 0)
        cat("\nNo misclassified quadrats\n\n")
    else {
        cat("\n", nmiss,
            "misclassified quadrats:\n\n")
        out <- data.frame("ind"=x$index, "Quadrat" = x$labels,
                          "Class" = x$class, "Predicted" = x$predicted,
                          "Division" = x$division)
        names(out)[1] <- ""
        print(out, row.names = FALSE)
    }
    invisible(x)
}
