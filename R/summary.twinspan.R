### summarize classification similarly as in the printed output of
### TWINSPAN (but in more compact form and only the essential
### information). The function traverses the classification tree
### recursively (depth-first).

#' Summary of Twinspan Classification
#'
#' The function gives a compact summary of divisions with indicator
#' species and items in final classification. The output gives the
#' same essential information as the printed output of TWINSPAN batch
#' program, but in more compact form.
#'
#' For each division, \code{summary} prints the eigenvalue, and for
#' quadrat divisions the indicator pseudospecies with their signs,
#' followed by \code{<} and the lowest score for the \dQuote{positive}
#' group. If the indicator score is below this value, follow the list
#' to the next item at the lower level, and to the second alternative
#' at the same level if the indicator score is at the limit or
#' higher. For division \eqn{k}, the next items are \eqn{2k}
#' (\dQuote{negative} group) or \eqn{2k+1} (\dQuote{positive}
#' group). For terminal groups, the function gives the size of the
#' group and lists its elements (quadrats or species).
#'
#' @return The function returns nothing. It only prints the result
#'     object in a human-readable way.
#'
#' @param object \code{\link{twinspan}} result object.
#' @param what Summarize either quadrat or species classification.
#' @param \dots Other arguments (ignored).
#'
#' @export
`summary.twinspan` <-
    function(object, what = c("quadrat","species"), ...)
{
    what <- match.arg(what)
    obj <- object[[what]]
    clid <- cut(object, what=what)
    len <- length(obj$eig) * 2 + 1
    state <- character(len)
    state[which(obj$eig > 0)] <- "division"
    state[unique(clid)] <- "cluster"
    ## twinvisit is called recursively
    o <- twinvisit(1, state, obj)
}

`twinvisit` <-
    function(k, state, obj)
{
    if(k > length(state) || state[k]=="")
        return(NULL)
    twinreport(k, state, obj)
    twinvisit(2*k, state, obj)
    twinvisit(2*k+1, state, obj)
}

`twinreport` <-
    function(k, state, obj)
{
    level <- sum(2^(1:10) <= k)
    cat(rep("  ", level), sep="")
    cat(k, ") ", sep="")
    if (state[k] == "division") {
        cat("eig=", round(obj$eig[k], 3), sep = "")
        if (!is.null(obj$indicators)) {
            ind <- obj$indicators[,k]
            ind <- ind[ind != 0]
            nm <- paste0(c("-","+")[(ind > 0) + 1], obj$indlabels[abs(ind)])
            cat(": ", nm)
            cat(" <", obj$positivelimit[k])
        }
        cat("\n")
    }
    else { # class
        nm <- obj$labels[obj$iclass == k]
        cat("N=", length(nm), ": ", sep="")
        cat(nm, "\n")
    }
}
