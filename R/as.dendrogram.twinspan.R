#' Extract Species or Quadrat Dendrograms
#'
#' Function extracts the species or quadrat classification as a
#' hierarchic \code{\link[stats]{dendrogram}}.
#'
#' The dendrogram heights are levels of divisions, total Chi-squares
#' of divisions and groups, or eigenvalues of divisions depending on
#' argument \code{height}. Terminal groups have no eigenvalues,
#' because they were not considered for division. For them the method
#' uses arbitrary value that for a group of \eqn{n} units is
#' proportion \eqn{(n-1)/n} of the height of mother
#' division. Chi-squares are evaluated also for terminal groups.
#' There is no guarantee that eigenvalues or Chi-squares decrease in
#' divisions, and there may be reversals where lower levels are higher
#' than their mother groups, and the plotted trees can be messy and
#' unreadable. Chi-squares decrease more monotonically than
#' eigenvalues of first axis.
#'
#' \R{} has a wealth of functions to handle and display
#' dendrograms. See \code{\link[stats]{dendrogram}} for general
#' description. There is even stronger support in packages (for
#' instance, \CRANpkg{dendextend}).
#'
#' The terminal groups of \code{\link{twinspan}} trees are not binary,
#' but may have several elements (quadrats, species).  In
#' \code{\link[stats]{dendrogram}} plots, it is best to set
#' \code{type="triangle"} for nicer looking trees.
#'
#' @return A \code{\link[stats]{dendrogram}} object.
#'
#' @examples
#'
#' ## Large datasets are difficult to show in dendrograms: take only
#' ## Northen Boreal quadrats (from 1 to 87).
#'
#' data(ahti)
#' tw <- twinspan(ahti[1:87,])
#' den <- as.dendrogram(tw)
#' str(den, max.level = 4)
#' plot(den, type = "triangle", nodePar = list(lab.cex=0.6, pch=NA))
#' den <- as.dendrogram(tw, height="chi")
#' plot(den, type = "triangle", nodePar = list(lab.cex=0.6, pch=NA))
#'
#' @seealso \code{\link{as.hclust.twinspan}}, \code{\link{dendrogram}}.
#'
#' @param object \code{\link{twinspan}} result object.
#' @param what Return either a \code{"quadrat"} or \code{"species"}
#'     dendrogram.
#' @param height Use either division levels (\code{"level"}), total
#'     Chi-squares (\code{"chi"}) or eigenvalues of first axis
#'     (\code{"eigen"}) of division as dendrogram heights.
#'
#' @param \dots Other parameters to functions.
#'
#' @importFrom stats as.dendrogram
#'
#' @export
`as.dendrogram.twinspan` <-
    function(object, what = c("quadrat", "species"),
             height = c("level", "chi", "eigen"), ...)
{
    what <- match.arg(what)
    height <- match.arg(height)
    obj <- object[[what]]
    clid <- cut(object, what=what)
    len <- length(obj$eig) * 2 + 1
    state <- character(len)
    state[which(obj$eig > 0)] <- "branch"
    state[unique(clid)] <- "leaf"
    ## eigen: divisions have eigenvalue, for terminal group of size n
    ## use proportion (n-1)/n of mother division
    if (height == "eigen") {
        eig <- numeric(len)
        eig[seq_along(obj$eig)] <- obj$eig
        for(k in which(state == "leaf")) {
            nk <- sum(clid == k)
            eig[k] <- (1 - 1/nk) * eig[k %/% 2]
        }
    } else if (height == "chi") {
        eig <- twintotalchi(object, what = what)
    } else {
        pow2 <- 2^(0:(object$levelmax+1))
        hmax <- sum(max(which(nchar(state) >0 )) >= pow2) + 1
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
            if (height %in% c("chi","eigen")) {
                hi <- 0
            } else {
                hi <- hmax - sum(k >= pow2) - 1
            }
            for (i in seq_len(length(zk))) {
                attr(zk[[i]], "label") <- labs[i]
                attr(zk[[i]], "members") <- 1L
                attr(zk[[i]], "height") <- hi
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
        ## Divisions have eigenvalue, but ev is never evaluated for
        ## terminal groups ("leaf"). We use an arbitrary value: for
        ## group of size n proportion (n-1)/n of the eigenvalue of
        ## mother division.
        if (height %in% c("chi","eigen")) {
            attr(zk, "height") <- eig[k]
        }
        else
            attr(zk, "height") <- hmax - sum(k >= pow2)
        z[[as.character(k)]] <- zk
    }
    structure(z[[as.character(k)]], class="dendrogram")
}
