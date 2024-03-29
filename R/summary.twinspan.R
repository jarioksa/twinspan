### summarize classification similarly as in the printed output of
### TWINSPAN (but in more compact form and only the essential
### information). The function traverses the classification tree
### recursively (depth-first).

#' Summary of twinspan Classification
#'
#' The function gives a compact summary of divisions with indicator
#' species and items in final classification. The output gives the
#' same essential information as the printed output of TWINSPAN batch
#' program, but in more compact form.
#'
#' For each division, \code{summary} prints the eigenvalue. For
#' quadrat divisions, it also prints the indicator pseudospecies with
#' their signs, followed by \code{<} and the lowest indicator score
#' for the \sQuote{positive} (right) group. If the indicator score is
#' below this value, follow the summary to the next item at the lower
#' level, and if the indicator score is at the limit or higher, follow
#' to the second alternative.  For division number \eqn{k}, the next
#' items are either \eqn{2k}{2*k} (\sQuote{negative} group) or
#' \eqn{2k+1}{2*k+1} (\sQuote{positive} group). Function
#' \code{\link{plot.twinspan}} displays the division numbers in a
#' classification tree.
#'
#' For terminal groups, the function gives the size of the group and
#' lists its elements (quadrats or species).
#'
#' @seealso \code{\link{plot.twinspan}} displays the same structure
#'     visually.  Function \code{\link{predict.twinspan}} follows the
#'     summary strcture to predict the classification with indicator
#'     pseudospecies.
#'
#' @examples
#' data(ahti)
#' tw <- twinspan(ahti)
#' summary(tw)
#' summary(tw, "species")
#'
#' @return The function returns nothing. It only prints the result
#'     object in a human-readable way.
#'
#' @param object \code{\link{twinspan}} result object.
#' @param what Summarize either quadrat or species classification.
#' @param binname Use binary labels for divisions instead of decimal numbers.
#' @param \dots Other arguments (ignored).
#'
#' @export
`summary.twinspan` <-
    function(object, what = c("quadrat","species"), binname = FALSE, ...)
{
    what <- match.arg(what)
    obj <- object[[what]]
    clid <- cut(object, what=what)
    len <- length(obj$eig) * 2 + 1
    state <- character(len)
    state[which(obj$eig > 0)] <- "division"
    state[unique(clid)] <- "cluster"
    ## twinvisit is called recursively
    o <- twinvisit(1, state, obj, binname = binname)
}

`twinvisit` <-
    function(k, state, obj, binname)
{
    if(k > length(state) || state[k]=="")
        return(NULL)
    twinreport(k, state, obj, binname = binname)
    twinvisit(2*k, state, obj, binname = binname)
    twinvisit(2*k+1, state, obj, binname = binname)
}

`twinreport` <-
    function(k, state, obj, binname = binname)
{
    level <- sum(2^(seq_len(15)) <= k)
    cat(rep("  ", level), sep="")
    if (binname)
        cat(class2bin(k), ") ", sep="")
    else
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
