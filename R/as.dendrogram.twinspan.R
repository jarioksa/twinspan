`as.dendrogram.twinspan` <-
    function(object, what = c("quadrat", "species"), eigenheight = FALSE, ...)
{
    what <- match.arg(what)
    obj <- object[[what]]
    clid <- cut(object, what=what)
    len <- length(obj$eig) * 2
    state <- character(len)
    state[which(obj$eig > 0)] <- "branch"
    state[unique(clid)] <- "leaf"
    if (eigenheight) {
        eig <- obj$eig
        hmax <- 1
    } else {
        hmax <- sum(max(which(nchar(state) >0 )) >= 2^(0:10)) + 1
    }
    z <- list()
    for(k in rev(seq_along(state))) {
        if(nchar(state[k]) == 0)
            next
        if(state[k] == "leaf") {
            zk <- as.list(which(clid == k))
            attr(zk, "members") <- length(zk)
            attr(zk, "midpoint") <- (length(zk)-1)/2
            labs <- obj$labels[clid == k]
            if (eigenheight) {
                height <- 0
            } else {
                height <- hmax - sum(k >= 2^(0:10)) - 1
            }
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
        if (eigenheight) {
            attr(zk, "height") <- if(k==1) 1 else eig[k %/% 2]
        }
        else
            attr(zk, "height") <- hmax - sum(k >= 2^(0:10))
        z[[as.character(k)]] <- zk
    }
    structure(z[[as.character(k)]], class="dendrogram")
}
