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

#' @export
`print.twinord` <- function(x, ...)
{
    cat("\nCorrespondence Analysis of twinspan division", x$division, "\n")
    cat("using vegan::decorana with ira=1 and twinspan downweighting\n")
    NextMethod("print")
}
