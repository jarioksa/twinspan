#' Summary of Twinspan Classification
#'
#' @export
`summary.twinspan` <-
    function(object, what = c("quadrat","species"), ...)
{
    what <- match.arg(what)
    obj <- object[[what]]
    clid <- cut(object, what=what)
    len <- length(obj$eig) * 4 + 1
    state <- character(len)
    state[which(obj$eig > 0)] <- "division"
    state[unique(clid)] <- "cluster"
    o <- twinvisit(1, state, obj)
}


`twinvisit` <-
    function(k, state, obj)
{
    if(state[k]=="")
        return(NULL)
    twinreport(k, state, obj)
    twinvisit(2*k, state, obj)
    twinvisit(2*k+1, state, obj)
}

`twinreport` <-
    function(k, state, obj)
{
    cat(k, ") ", sep="")
    if (state[k] == "division") {
        cat("eig=", round(obj$eig[k], 3), sep = "")
        if (!is.null(obj$indicators)) {
            ind <- obj$indicators[,k]
            ind <- ind[ind != 0]
            nm <- paste0(c("-","+")[ind > 0], obj$indlabels[abs(ind)])
            cat(": if ", nm)
            cat(" <", obj$positivelimit[k], "goto", 2*k, "else", 2*k+1)
        }
        cat("\n")
    }
    else { # class
        nm <- obj$labels[obj$iclass == k]
        cat("N=", length(nm), ":", nm, "\n")
    }
}
