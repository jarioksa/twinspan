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
#' returns only the index of misclassified quadrats. In default, it
#' also returns the names, twinspan classes and predicted classes of
#' misclassified quadrats, and the division where the
#' misclassification occurred and classifications diverged. The
#' divisions and their numbers can be seen in
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
#' @seealso \code{\link{cut.twinspan}},
#'     \code{\link{predict.twinspan}}, \code{\link{summary.twinspan}},
#'     \code{\link{plot.twinspan}}.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' misclassified(tw)
#'
#' @param x \code{\link{twinspan}} result object.
#' @param verbose If \code{FALSE}, returns only the index of
#'     misclassified quadrats.
#'
#' @export
`misclassified` <-
    function(x, verbose = TRUE)
{
    class <- cut(x)
    pred <- predict(x)
    # misclassified cases
    ind <- which(class != pred)
    if (!verbose)
        return(ind)
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
            "misclassified quadrats and their level of divergence:\n\n")
        out <- data.frame("ind"=x$index, "Quadrat" = x$labels,
                          "Class" = x$class, "Predicted" = x$predicted,
                          "Division" = x$division)
        names(out)[1] <- ""
        print(out, row.names = FALSE)
    }
    invisible(x)
}
