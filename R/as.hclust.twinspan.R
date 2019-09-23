### Quadrat divisions without terminal items (quadrats) as an hclust
### tree. Hclust trees can not have polytomies, and quadrats must be
### omitted. There already is as.dendrogram, but try to see if this
### can be developed into a decision tree.

#' Extract Quadrat Grouping as Hierarchical Cluster Tree
#'
#' Function extracts quadrat classification as an \code{\link{hclust}}
#' object. The terminal items are the quadrat classes, but quadrats
#' are not shown: \code{\link{hclust}} cannot handle polytomies that
#' are needed to display individual quadrats.  Use
#' \code{\link{as.dendrogram}} to show the single items or species
#' classification.
#'
#' @param x \code{\link{twinspan}} result object.
#' @param \dots Other parameters to the function (ignored).
#'
#' @importFrom stats as.hclust
#'
#' @export
`as.hclust.twinspan` <-
    function(x, ...)
{
    class <- cut(x)
    n <- max(class)
    state <- character(n)
    state[which(x$quadrat$eig > 0)] <- "div"
    state[unique(class)] <- "clust"
    nclus <- sum(state == "clust")
    line <- numeric(n)
    merge <- matrix(0, nclus-1, 2)
    height <- numeric(nclus-1)
    pow2 <- 2^(1:15)
    maxh <- sum(max(class) >= pow2) + 1
    nro <- nclus
    now <- 0
    for(i in seq(n-1, 2, by=-2)) {
        if (state[i] == "")
            next
        now <- now + 1
        height[now] <- maxh - sum(i >= pow2)
        for(j in 2:1) {
            line[i+j-1] <- now
            if(state[i+j-1] == "clust") {
                merge[now,j] <- -nro
                nro <-  nro - 1
            } else {
                merge[now,j] <- line[(i + j - 1) * 2]
            }
        }
    }
    labels <- table(class)
    labels <- paste0(names(labels), " (N=", labels, ")")
    nodelabels <- rev(which(state=="div"))
    ind <- x$quadrat$index
    order <- order(tapply(order(ind), class, min))
    out <- list(merge = merge, labels = labels, height = height, order = order,
                nodelabels = nodelabels, method = "twinspan")
    class(out) <- "hclust"
    out
}

### plot.twinspan as plot of hclust tree

#' Plot Classificaton Tree of Quadrats
#'
#' @param x \code{\link{twinspan}} result object.
#' @param main Main title of the plot.
#' @param \dots Other parameters passed to \code{\link{plot}} and
#'     \code{\link[vegan]{ordilabel}}.
#'
#' @importFrom vegan ordilabel
#' @importFrom graphics plot
#'
#' @export
`plot.twinspan` <-
    function(x, main = "Twinspan Dendrogram", ...)
{
    x <- as.hclust(x)
    plot(x, main = main, ...)
    ordilabel(x, "internal", labels = x$nodelabels, ...)
}
