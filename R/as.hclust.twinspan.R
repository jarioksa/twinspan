### Quadrat divisions without terminal items (quadrats) as an hclust
### tree. Hclust trees can not have polytomies, and quadrats must be
### omitted. There already is as.dendrogram, but try to see if this
### can be developed into a decision tree.

#' Extract twinspan Grouping as Hierarchical Cluster Tree
#'
#' Function extracts classification as an \code{\link{hclust}}
#' object. The terminal items are the final groups, but quadrats or
#' species are not shown: \code{\link{hclust}} cannot handle
#' polytomies that are needed to display group members.  Use
#' \code{\link{as.dendrogram}} to show the single items. The group ID
#' number and number of items in the terminal group are used as group
#' names and are displayed in plots.
#'
#' @seealso \code{\link{as.dendrogram.twinspan}},
#'     \code{\link{hclust}}, \code{\link{plot.twinspan}},
#'     \code{\link{image.twinspan}}.
#'
#' @return an \code{\link{hclust}} object.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' plot(as.hclust(tw))
#' plot(as.hclust(tw, "species"))
#'
#' @param x \code{\link{twinspan}} result object.
#' @param what Extract \code{"quadrat"} or \code{"species"}
#'     classification tree.
#' @param height Use either division levels (\code{"level"}) or total
#'     Chi-squares of division (\code{"chi"}) as heights of internal
#'     nodes in the tree.
#' @param \dots Other parameters to the function (ignored).
#'
#' @importFrom stats as.hclust
#'
#' @export
`as.hclust.twinspan` <-
    function(x, what = c("quadrat","species"), height = c("level", "chi"),
             ...)
{
    what <- match.arg(what)
    height <- match.arg(height)
    class <- cut(x, what=what)
    n <- max(class)
    state <- character(n)
    state[which(x[[what]]$eig > 0)] <- "div"
    state[unique(class)] <- "clust"
    nclus <- sum(state == "clust")
    line <- numeric(n)
    merge <- matrix(0, nclus-1, 2)
    treeheight <- numeric(nclus-1)
    pow2 <- 2^(seq_len(x$levelmax+1))
    maxh <- sum(max(class) >= pow2) + 1
    nro <- nclus
    now <- 0
    for(i in seq(n-1, 2, by=-2)) {
        if (state[i] == "")
            next
        now <- now + 1
        treeheight[now] <- maxh - sum(i >= pow2)
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
    if (height == "chi") {
        ## replace levels with Chi-squares of internal nodes
        treeheight <- twintotalchi(x, what)[nodelabels]
        ## height should be ordered in hclust tree so that hclust
        ## methods know how to handle tree
        if (is.unsorted(treeheight)) {
            o <- order(treeheight)
            oo <- order(o)
            for (i in 1:nrow(merge))
                for (j in 1:2)
                    if (merge[i,j] > 0) # a (reordered) division
                        merge[i,j] <- oo[merge[i,j]]
            merge <- merge[o,]
            treeheight <- treeheight[o]
            nodelabels <- nodelabels[o]
        }
    }
    ind <- x[[what]]$index
    order <- order(tapply(order(ind), class, min))
    out <- list(merge = merge, labels = labels, height = treeheight,
                order = order, nodelabels = nodelabels, method = "twinspan")
    class(out) <- "hclust"
    out
}

### plot.twinspan as plot of hclust tree

#' Plot Classification Tree
#'
#' Function displays the classification tree.
#'
#' The internal nodes are labelled by the numbers of division. These
#' are the same numbers as used in \code{\link{summary.twinspan}} and
#' returned by \code{\link{cut.twinspan}} or
#' \code{\link{predict.twinspan}} for the same classification
#' level. For terminal groups the plot shows the number of the group
#' and the number of items (quadrats or species) in the group. For
#' division number \eqn{k}, its daughter divisions or groups are
#' \eqn{2k}{2*k} and \eqn{2k+1}{2*k+1}. The tree is similar as a plot
#' of \code{\link{as.hclust.twinspan}}, but adds numbers of internal
#' nodes.
#'
#' @seealso \code{\link{summary.twinspan}} for similar textual
#'     presentation also showing the items (quadrats, species) in
#'     terminal groups. \CRANpkg{vegan} function
#'     \code{\link[vegan]{scores.hclust}} can extract the coordinates
#'     of internal (or terminal nodes), and
#'     \code{\link[vegan]{ordilabel}} is used add the labels on
#'     internal nodes.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' plot(tw, "species")
#'
#' @param x \code{\link{twinspan}} result object.
#' @param what Plot \code{"quadrat"} or \code{"species"}
#'     classification tree.
#' @param main Main title of the plot.
#' @param \dots Other parameters passed to \code{\link{plot}} and
#'     \code{\link[vegan]{ordilabel}}.
#' @param height Use either division levels (\code{"level"}) or total
#'     Chi-squares of division (\code{"chi"}) as heights of internal
#'     nodes in the tree.
#' @importFrom vegan ordilabel
#' @importFrom graphics plot
#'
#' @export
`plot.twinspan` <-
    function(x, what = c("quadrat", "species"), height = c("level", "chi"),
             main = "Twinspan Dendrogram",
             ...)
{
    what <- match.arg(what)
    x <- as.hclust(x, what = what, height = height)
    plot(x, main = main, ...)
    ordilabel(x, "internal", labels = x$nodelabels, ...)
}
