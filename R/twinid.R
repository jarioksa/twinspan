##' Get Twinspan Classs Identiers for Clustering Object
##'
##' @description Twinspan returns the classification topology as a
##'     single integer vector. These functions find similar
##'     classification identifiers for each sampling unit, or
##'     \code{cut} that vector for a lower number of classes.
##'
##' @details Twinspan expresses the topology of cluster tree as an
##'     integer. When a cluster \eqn{z} is split into two, its
##'     daughters will be \eqn{2z}{2*z} and \eqn{2z+1}{2*z+1}, and its
##'     parent cluster is found with integer division \eqn{z/2}. The
##'     classification vector only stores the topology of the trees,
##'     and has no information on heights.
##'
##'     \code{\link{twinspan}} will not split small clusters and only
##'     proceeds to a defined depth of divisions. In contrast,
##'     \code{twinind} proceeds to each terminal unit (leaf, sampling
##'     unit, quadrat) and these will all have unique
##'     identifiers. With \code{cut} function you can restrict the
##'     identifiers to certain level of classification similarly as in
##'     \code{twinspan} (see \code{\link{cut.twinspan}}).
##'
##' @section Warning: If the classification is deep and has many (>
##'     30) levels of hierarchy, the identifiers can exceed the
##'     integer maximum in \R{}, and leaves may have non-unique
##'     identifiers, and may not recover the correct
##'     topoloty. However, they may still be unique beyond this limit,
##'     but the user should check this after getting a warning.
##'
##' @return Vector of class \code{"twinid"} giving
##'     \code{\link{twinspan}} id of each sampling unit.
##'
##' @param hclust Cluster Analysis result compatible with
##'     \code{\link{hclust}}.
##' @param x Vector of classification IDs from \code{twinind}.
##' @param level Level of hierarchy of classification. If missing,
##'     level used in the object will be returned.
##' @param ... Other parameters to functions (ignored).
##'
##' @examples
##' data(ahti)
##' cl <- hclust(dist(ahti, "manhattan"), "average")
##' (id <- twinid(cl))
##' cut(id, 6)
##' table(cut(id, 6))
##'
##' @export
`twinid` <-
    function(hclust, ...)
{
    ## Works with hclust clusters. These have a matrix of dims (n-1)
    ## times 2, where each row defines a fusion of two items. Negative
    ## entries -1...-n denote original observations (leaves, terminal
    ## nodes) and positive entries previous rows of the merge matrix
    ## (internal nodes). Twinspan in contrast defines topology by an
    ## integer vector, where the daughters of a mother cluster x are
    ## 2*x and 2*x+1, and the mother cluster of id z can be found by
    ## integer division z %/% 2 which allows recontruction of the tree
    ## topology from a single integer vector.
    if (!inherits(hclust, "hclust"))
        stop("'hclust' must be a clustering object from hclust() function")
    merge <- hclust$merge
    ## id, leftid, classid are "global" variables which will get
    ## superassigned (<<-) values in visit().
    id <- 1
    leftid <- numeric(nrow(merge))
    classid <- numeric(nrow(merge) + 1)
    ## Recursive function to label units (leaves) as Twinspan class IDs
    visit <- function(i, j) {
        if (j==1) {
            ## left: deeper level, double id and save on left identifiers
            id <<- 2*id
            leftid[i] <<- id
        } else {
            ## right: same level as left, but increase id by one
            id <<- leftid[i] + 1
        }
        if (merge[i,j] < 0) {
            ## leaf: save its id and exit this instance of visit()
            classid[-merge[i,j]] <<- id
        } else {
            ## internal node: visit recursively left and right branches
            visit(merge[i, j], 1)
            visit(merge[i, j], 2)
        }
    }
    visit(nrow(merge), 1)
    visit(nrow(merge), 2)
    ## all leaves should have unique identifiers
    if (any(duplicated(classid)))
        stop("some terminal nodes (leaves) have duplicated identifiers")
    ## warn on large IDs: depth > 30 can exceed integer maximum, and
    ## we cannot guarantee integer division (even with unique IDs)
    if (max(classid) > .Machine$integer.max)
        warning("some class identifiers larger than integer maximum")
    class(classid) <- "twinid"
    classid
}

#' @rdname twinid
#' @export
`cut.twinid` <-
    function(x, level, ...)
{
    ## no level: return x
    if (missing(level))
        return(x)
    ## max id for given level
    clmax <- 2^(level+1) - 1
    while (any(big <- x > clmax)) {
        ## mother class by integer division
        x[big] <- x[big] %/% 2
    }
    x
}
