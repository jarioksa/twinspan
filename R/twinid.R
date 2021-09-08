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
    classid
}
