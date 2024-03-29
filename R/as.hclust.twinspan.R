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
#' Function can return either a tree showing the \code{twinspan}
#' hierarchy or showing the heterogeneity of each group or
#' division. In the first case, all divisions and groups at a certain
#' level of hierarchy are at the same height, but in the latter the
#' divisions are at the height defined by their heterogeneity. The
#' criterion of heterogeneity is the total chi-square (also known as
#' inertia) of the matrix that \code{\link{twinspan}} internally uses
#' in that division (see \code{\link{twintotalchi}}). This tree gives
#' the visual presentation of the modified method of Roleček et
#' al. (2009).
#'
#' When tree heights are based on heterogeneity, subgroups can be more
#' heterogeneous than their parent group. These appear as reversed
#' branches in the tree. A warning is issued for each such case.
#'
#' @encoding UTF-8
#'
#' @seealso \code{\link{as.dendrogram.twinspan}} provides an
#'     alternative which also shows the sampling units (quadrats or
#'     species). The result is based \code{\link{hclust}} and can be
#'     handled with its support
#'     functions. \code{\link{plot.twinspan}},
#'     \code{\link{image.twinspan}} display the tree. Function
#'     \code{\link{cut.twinspan}} cuts the tree by a level of
#'     hierarchy, and \code{\link{cuth}} by heterogeneity for original
#'     sampling units (quadrats, species).
#'
#' @return an \code{\link{hclust}} object amended with labels for
#'     internal nodes (\code{nodelabels}).
#'
#' @references
#' Roleček, J, Tichý, L., Zelený, D. & Chytrý, M. (2009). Modified
#' TWINSPAN classification in which the hierarchy respects cluster
#' heterogeneity. \emph{J Veg Sci} 20: 596--602.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' plot(as.hclust(tw, "species"))
#' cl <- as.hclust(tw)
#' ## plot and 8 groups by hierarchy level
#' plot(cl)
#' rect.hclust(cl, 8)
#' ## plot and 8 groups by heterogeneity
#' cl <- as.hclust(tw, height="chi")
#' plot(cl)
#' rect.hclust(cl, 8)
#'
#'
#'
#' @param x \code{\link{twinspan}} result object.
#' @param what Extract \code{"quadrat"} or \code{"species"}
#'     classification tree.
#' @param height Use either division levels (\code{"level"}) or total
#'     Chi-squares of division (\code{"chi"}) as heights of internal
#'     nodes in the tree.
#' @param binname Use binary labels for classes instead of decimal
#'     numbers.
#' @param \dots Other parameters to the function (ignored).
#'
#' @importFrom stats as.hclust
#'
#' @export
`as.hclust.twinspan` <-
    function(x, what = c("quadrat","species"), height = c("level", "chi"),
             binname = FALSE, ...)
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
    if (binname)
        names(labels) <- sapply(as.numeric(names(labels)), class2bin)
    labels <- paste0(names(labels), " (N=", labels, ")")
    nodelabels <- rev(which(state=="div"))
    if (height == "chi") {
        ## replace levels with Chi-squares of internal nodes
        treeheight <- twintotalchi(x, what)[nodelabels]
        ## height should be ordered in hclust tree so that hclust
        ## methods know how to handle tree
        if (is.unsorted(treeheight)) {
            o <- order(treeheight)
            o <- o[rev(fixTreeReversal(rev(nodelabels[o]), index = TRUE))]
            oo <- order(o)
            for (i in 1:nrow(merge))
                for (j in 1:2)
                    if (merge[i,j] > 0) # a (reordered) division
                        merge[i,j] <- oo[merge[i,j]]
            merge <- merge[o,]
            treeheight <- treeheight[o]
            nodelabels <- nodelabels[o]
            ## tree reversals should be handled, but check anyway
            if (any(merge > row(merge))) # mother before her daughters
                stop("some split groups are more heterogeneous than their mother")
        }
    }
    ind <- x[[what]]$index
    order <- order(tapply(order(ind), class, min))
    out <- list(merge = merge, labels = labels, height = treeheight,
                order = order, nodelabels = nodelabels, method = "twinspan")
    class(out) <- "hclust"
    out
}

## If within-group heterogeneity (such as Chi-square) is used as
## height in a twinspan tree, there may be reversed branches or the
## subgroup is more heterogeneous than her mother group. This function
## goes through order of groups and if it finds a subgroup was used
## before the parent group, it will move the kid behind her
## parent. Non-exported support function tuned to work only with our
## as.hclust.twinspan and cuth functions.

fixTreeReversal <-
    function(order, index = FALSE)
{
    n <- length(order)
    ## parent of group i is i %/% 2
    if(index)
        idx <- rev(seq_along(order)) # reverse index
    repeat{
        ## may need several passages
        SWAPPED <- FALSE
        for (i in 2:n) {
            a <- order[1:(i-1)] %/% 2 == order[i]
            if (any(a)) {
                SWAPPED <- TRUE
                k <- min(which(a)) # there should be no two kids, but check
                ## move kid immediately after her parent
                if (index) {
                    idx[k:i] <- idx[c((k+1):i, k)]
                }
                order[k:i] <- order[c((k+1):i, k)]
                warning(
                    gettextf("tree reversal: group %d more heterogeneous than parent %d",
                             order[i], order[i-1]),
                    call. = FALSE)
            }
        }
        if (!SWAPPED) break
    }
    if(index) idx else order
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
#' level. For terminal groups the plot shows the numeric code of the
#' group and the number of items (quadrats or species) in the
#' group. For division number \eqn{k}, its daughter divisions or
#' groups are coded \eqn{2k}{2*k} and \eqn{2k+1}{2*k+1}. The tree is
#' similar as a plot of \code{\link{as.hclust.twinspan}}, but adds
#' numbers of internal nodes. The tree can be based either on the
#' levels of hierarchy or on the heterogeneity of division as assessed
#' by chi-square (or inertia) of the division (see
#' \code{\link{as.hclust.twinspan}}).
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
#' ## default plot for quadrats
#' plot(tw)
#' ## plot by the heterogeneity of divisions
#' plot(tw, height = "chi")
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
#' @param binname Use binary labels for classes and nodes instead of
#'     decimal numbers.
#' @importFrom vegan ordilabel
#' @importFrom graphics plot
#'
#' @export
`plot.twinspan` <-
    function(x, what = c("quadrat", "species"), height = c("level", "chi"),
             main = "Twinspan Dendrogram", binname = FALSE, ...)
{
    what <- match.arg(what)
    x <- as.hclust(x, what = what, height = height, binname = binname)
    plot(x, main = main, ...)
    labels <- x$nodelabels
    if (binname)
        labels <- sapply(labels, class2bin)
    ordilabel(x, "internal", labels = labels, ...)
}
