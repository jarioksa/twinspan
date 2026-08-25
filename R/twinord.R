### Ordination of twinspan division

#' Ordination Demonstrating Twinspan Division
#'
#' Function performs orthogonal correspondence analysis with
#' downweighting for quadrats in one division (\pkg{vegan} function
#' \code{\link[vegan]{decorana}}). This analysis is similar as the one
#' performed in \code{twinspan} when it splits the division. The
#' results can be used to visualize and analyse the division. Function
#' returns the ordination object with auxiliary information on
#' subdivision, quadrat indicator scores, indicator pseudospecies and
#' indicator values of pseudospecies.
#'
#' @return Function returns an object of class \code{"twinord"} which
#'     inherits from a \pkg{vegan} \code{\link[vegan]{decorana}}
#'     object, but adds the following auxiliary elements:
#'
#' \describe{
#' \item{class}{Two classes into which this division was split.}
#' \item{indscores}{Indicator scores for quadrats from \code{\link{indscores}}.}
#' \item{indicators}{Name labels of indicator pseudospecies used to split this
#'   class from \code{\link{twinspan}} (see \code{\link{summary.twinspan}}).}
#' \item{indvals}{Indicator values for pseudospecies from
#'   \code{\link{indvalues}}.}
#' \item{division}{\code{twinspan} division.}
#' }
#'
#' @param object \code{twinspan} result object.
#' @param division Number code of \code{twinspan} division for quadrats.
#'
#' @importFrom stats cor
#' @importFrom vegan decorana
#' @export
`twinord` <-
    function(object, division)
{
    if(object$quadrat$eig[division] == 0 ||
       length(object$quadrat$eig) < division)
        stop(division, " is not a division")
    ## ordination
    ord <- decorana(twin2stack(object, subset = twingroup(object, division),
                               downweight = TRUE),
                    ira = 1, iresc = 0)
    ## subdivisions
    level <- floor(log2(division))
    cl <- cut(object, level + 1)[twingroup(object, division)]
    ## indicator scores for quadrats
    indscore <- indscores(object, division)$score
    ## take care that negative group is on negative side of axis 1
    if (cor(indscore, ord$rproj[,1]) < 0) {
        ord$rproj[,1] <- -ord$rproj[,1]
        ord$cproj[,1] <- -ord$cproj[,1]
    }
    ## indicator species
    inds <- object$quadrat$indlabels[abs(object$quadrat$indicators[,division])]
    ## indicator values for species
    indval <- indvalues(object, division)
    ## add auxiliary items to the decorana object
    ord$class <- cl
    ord$indscores <- indscore
    ord$indicators <- inds
    ord$indvalues <- indval
    ord$division <- division
    ord$call <- match.call()
    class(ord) <- c("twinord", class(ord))
    ord
}

### non-document print adds two lines to print.decorana
#' @export
`print.twinord` <- function(x, ...)
{
    cat("\nCorrespondence Analysis of twinspan division", x$division, "\n")
    cat("using vegan::decorana with ira=1 and twinspan downweighting\n")
    NextMethod("print")
}

### plot

#' @rdname twinord
#'
#' @param x \code{twinord} result object.
#' @param quadrats Plot quadrats as text, points or with their
#'     indicator score or none.
#' @param pseudospecies Plot text for none, all, indicator
#'     pseudospecies or valuable species exceeding \code{value}.
#' @param value Limit of absolute indicator value for displayed
#'     pseudospecies.
#' @param qcol Colour of quadrat text in negative and positive groups.
#' @param \dots Other arguments passed to graphical functions.
#'
#' @importFrom graphics points text
#' @importFrom vegan ordihull
#'
#' @export
`plot.twinord` <- function(x, quadrats = c("text", "points", "score", "none"),
                           pseudospecies =
                               c("none", "all", "indicators", "valuable"),
                           value = 0.25,
                           qcol = c("blue", "red"), ...)
{
    quadrats <- match.arg(quadrats)
    pseudospecies <- match.arg(pseudospecies)
    cl12 <- x$class - 2 * x$division + 1
    ## Handle species data. For appropriate scaling it is best adjust
    ## data in input x
    if (pseudospecies == "indicators")
        x$cproj <- x$cproj[x$indicators, , drop=FALSE]
    else if (pseudospecies == "valuable") {
        i <- rev(order(abs(x$indvalues)))
        indval <- x$indvalues[i]
        x$cproj <- x$cproj[i,]
        val <- abs(indval) >= value
        ## only the best pseudospecies of each species
        nodup <- !duplicated(gsub("[1-9]$", "", rownames(x$cproj)))
        x$cproj <- x$cproj[val & nodup, , drop=FALSE]
    }

    display <- if (pseudospecies == "none" && quadrats == "none")
                   "none"
               else if (pseudospecies == "none")
                   "sites"
               else if (quadrats == "none")
                   "species"
               else
                   "both"

    pl <- vegan:::plot.decorana(x, display = display, type = "n", ...)
    switch(quadrats,
           "text" = text(pl, what = "sites", col = qcol[cl12],
                         optimize = TRUE, xpd = TRUE, ...),
           "points" = points(pl, what = "sites", col = qcol[cl12], ...),
           "score" = text(pl, what = "sites", labels = x$indscores,
                          col = qcol[cl12], ...),
           "none" = NULL
           )
    if (quadrats != "none")
        ordihull(x, x$class, draw = "polygon", label = TRUE, cex=1)
    if (pseudospecies != "none")
        text(pl, what = "species", optimize=TRUE, xpd = TRUE, col = 1, ...)
    invisible(pl)
}
