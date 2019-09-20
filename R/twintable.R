### TWINSPAN table

#' Community Table Ordered by Twinspan Classification
#'
#' Prints a community table of pseudospecies ordered by
#' \code{\link{twinspan}} classification.
#'
#' @param object \code{\link{twinspan}} result object.
#' @param maxspp Maximum number of most abundant species
#'     displayed. The abundance is estimated with pseudospecies cut
#'     levels. The default is to show all species.
#'
#' @importFrom vegan vegemite
#'
#' @export
`twintable` <-
    function(object, maxspp)
{
    i <- object$quadrat$index
    j <- object$species$index
    mat <- twin2mat(object)
    ## take only maxspp most abundant species - but if there are ties,
    ## take all tied species, even if we go over maxspp
    if (!missing(maxspp) && ncol(mat) > maxspp) {
        sptot <- colSums(mat)
        maxlim <- sort(sptot, decreasing = TRUE)[maxspp]
        j <- j[sptot[j] >= maxlim]
    }

    vegemite(mat[i,j], zero="-")
}
