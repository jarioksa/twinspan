### Twinspan ID for clustering object
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
    ## warn on large IDs: depth > 30, floor(log2(classid)) gives the depths.
    if (max(classid) > .Machine$integer.max)
        warning("some class identifiers larger than integer maximum")
    class(classid) <- "twinid"
    classid
}

`cut.twinid` <-
    function(x, level, ...)
{
    ## max id for given level
    clmax <- 2^(level+1) - 1
    while (any(big <- x > clmax)) {
        ## mother class by integer division
        x[big] <- x[big] %/% 2
    }
    x
}
