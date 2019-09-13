`as.dendrogram.twinspan` <-
    function(object, ...)
{
    clid <- cut(object)
    len <- length(object$quadrat$eig)
    state <- character(len)
    state[unique(clid)] <- "leaf"
    state[object$quadrat$eig > 0] <- "branch"
    z <- list()
    for(k in rev(seq_len(len))) {
        if(nchar(state[k]) == 0)
            next
        if(state[k] == "leaf") {
            zk <- as.list(which(clid == k))
            attr(zk, "members") <- length(zk)
            attr(zk, "midpoint") <- (length(zk)-1)/2
            labs <- object$quadrat$labels[clid == k]
            height <- 6 - sum(k >= 2^(0:7))
            for (i in seq_len(length(zk))) {
                attr(zk[[i]], "label") <- labs[i]
                attr(zk[[i]], "members") <- 1L
                attr(zk[[i]], "height") <- height
                attr(zk[[i]], "leaf") <- TRUE
            }
        }
        else { # a branch
            x <- c(2*k, 2*k+1)
            x <- as.character(x)
            zk <- list(z[[x[1]]], z[[x[2]]])
            attr(zk, "members") <- attr(z[[x[1]]], "members") +
                attr(z[[x[2]]], "members")
            attr(zk, "midpoint") <- (attr(z[[x[1]]], "members") +
                                     attr(z[[x[1]]], "midpoint") +
                                     attr(z[[x[2]]], "midpoint"))/2
            z[[x[1]]] <- z[[x[2]]] <- NULL
        }
        attr(zk, "height") <- 7 - sum(k >= 2^(0:7))
        z[[as.character(k)]] <- zk
    }
    structure(z[[as.character(k)]], class="dendrogram")
}
