### tabasco function to draw a classification heatmap with species and
### quadrat trees on the margins, and cells showing the average
### abundance in species x quadrat cell. Function basically uses chain
### vegan::tabasco -> stats::heatmap -> graphics::image, and since
### image is a generic function, we can have method for it.

#' Image of Quadrat and Species Classification and Average Abundances in Cells
#'
#' Function draws an \code{\link{image}} of mean pseudospecies values
#' in cells of final quadrat and species classifications together with
#' corresponding cluster trees.
#'
#' The function is based on \CRANpkg{vegan} function
#' \code{\link[vegan]{tabasco}}.  Mean abundances in cells are shown
#' by colours which can be modified by the user (see
#' \code{\link[vegan]{tabasco}}). The species groups are rows and
#' quadrat groups are columns. These are labelled by their decimal
#' group numbers with information of group sizes. These numbers and
#' the order of the groups are similar as in hierarchic cluster trees
#' produced by \code{\link{as.hclust.twinspan}} and
#' \code{\link{summary.twinspan}}. Function also adds quadrat
#' classification tree on the top and species classification tree on
#' the left margin (see \code{\link{as.hclust.twinspan}}).
#'
#' @return Function returns invisibly the matrix of mean abundances
#'     used to produce the graph.
#'
#' @param x \code{\link{twinspan}} result object.
#' @param \dots Other arguments passed to \code{\link[vegan]{tabasco}}
#'     and further to \code{\link[stats]{heatmap}}.
#'
#' @importFrom graphics image
#' @importFrom stats as.hclust
#' @importFrom vegan tabasco
#'
#' @export
`image.twinspan` <-
    function(x, ...)
{
    mat <- twin2mat(x) # matrix of pseudospecies data
    spcl <- as.hclust(x, "species") # species tree
    qcl <- as.hclust(x) # quadrat tree
    ## calculate means of species cluster x quadrat cluster cells
    mat <- apply(mat, 1, function(z)
        tapply(z, cut(x, what="species"), mean))
    mat <- apply(mat, 1, function(z)
        tapply(z, cut(x, what="quadrat"), mean))
    rownames(mat) <- qcl$labels
    colnames(mat) <- spcl$labels
    tabasco(mat, use=qcl, sp.ind = spcl, Rowv = FALSE, Colv = FALSE, ...)
    invisible(mat)
}
