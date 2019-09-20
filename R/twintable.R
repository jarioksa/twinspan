### TWINSPAN table

#' Community Table Ordered by Twinspan Classification
#'
#' Prints a community table of pseudospecies ordered by
#' \code{\link{twinspan}} classification.
#'
#' @param object \code{\link{twinspan}} result object.
#'
#' @importFrom vegan vegemite
#'
#' @export
`twintable` <-
    function(object)
{
    i <- object$quadrat$index
    j <- object$species$index
    mat <- twin2mat(object)
    vegemite(mat[i,j], zero="-")
}
