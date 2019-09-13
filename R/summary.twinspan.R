#' Summary of Twinspan Classification
#'
#' @export
`summary.twinspan` <-
    function(object, what = c("quadrat", "species"), ...)
{
    what <- match.arg(what)
    obj <- object[[what]]
    clid <- cut(object, what=what)
    len <- length(obj$eig) * 2 + 1
    state <- character(len)
    state[which(obj$eig > 0)] <- "division"
    state[unique(clid)] <- "cluster"
    z <- list()
    for(k in rev(seq_along(state))) {
        if(nchar(state[k]) == 0)
            next
        if(state[k] == "cluster") {
            items <- obj$labels[clid==k]
            N <- length(items)
            zk <- list(type = "twinclus",
                       items = items,
                       N = N, k = k)
        }
        else { # a division
            kids <- c(2*k, 2*k+1)
            kids <- as.character(kids)
            N <- z[[kids[1]]]$N + z[[kids[2]]]$N
            eig <- obj$eig[k]
            i <- obj$indicators[,k]
            i <- i[i != 0]
            inds <- paste0(c("-","+")[(i>0)+1], obj$indlabels[abs(i)])
            positivelimit <- obj$positivelimit[k]
            zk <- list(type = "twindiv", N = N, eig = eig, indicator = inds,
                       positivelimit = positivelimit, k = k,
                       kids=list(z[[kids[1]]], z[[kids[2]]]))
            z[[kids[1]]] <- z[[kids[2]]] <- NULL
        }
        z[[as.character(k)]] <- zk
    }
    structure(z[[as.character(k)]], class="summary.twinspan")
}
