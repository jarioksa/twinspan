### Twinspan ID for clustering object
`twinid` <-
    function(x, ...)
{
    merge <- x$merge
    ## id, leftid, classid are "global" variables which will get
    ## superassigned (<<-) values in visit()
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
