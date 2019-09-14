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
    o <- twinvisit(1, state)
}


`twinvisit` <-
    function(k, state)
{
    if(state[k]=="")
        return(NULL)
    twinreport(k, state)
    twinvisit(2*k, state)
    twinvisit(2*k+1, state)
}

`twinreport` <-
    function(k, state)
{
    cat(k, ") ", state[k], "\n", sep="")
}
