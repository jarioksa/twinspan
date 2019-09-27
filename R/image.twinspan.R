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
#' \code{\link[vegan]{tabasco}}).  The mean abundances are either the
#' average values of all species of the species group (default), or
#' the average of the most abundant species of the species group (tied
#' abundances broken by species frequency) if
#' \code{leadingspecies=TRUE}. The species groups are rows and quadrat
#' groups are columns. These are labelled by their decimal group
#' numbers with information of group sizes, or species by their name
#' if only leading species are shown. These numbers and the order of
#' the groups are similar as in hierarchic cluster trees produced by
#' \code{\link{as.hclust.twinspan}} and
#' \code{\link{summary.twinspan}}.  Function also adds quadrat
#' classification tree on the top and species classification tree on
#' the left margin (see \code{\link{as.hclust.twinspan}}).
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' image(tw)
#' image(tw, leading = TRUE)
#'
#' @return Function returns invisibly the matrix of mean abundances
#'     used to produce the graph.
#'
#' @param x \code{\link{twinspan}} result object.
#' @param leadingspecies Show averages of leading species (most
#'     abundant, ties broken by frequency) of each cluster instead of
#'     averages of all species of the group.
#' @param \dots Other arguments passed to \code{\link[vegan]{tabasco}}
#'     and further to \code{\link[stats]{heatmap}}.
#'
#' @importFrom graphics image
#' @importFrom stats as.hclust
#' @importFrom vegan tabasco
#'
#' @export
`image.twinspan` <-
    function(x, leadingspecies = FALSE, ...)
{
    mat <- twin2mat(x) # matrix of pseudospecies data
    spcl <- as.hclust(x, "species") # species tree
    qcl <- as.hclust(x) # quadrat tree
    ## calculate quadrat group means of leading species
    if (leadingspecies) {
        j <- goodspec(x, "leading")
        spcl$labels <- colnames(mat)[j]
        mat <- mat[,j]
        mat <- apply(mat, 2, function(z)
            tapply(z, cut(x, what="quadrat"), mean))
    } else {
        ## calculate means of species cluster x quadrat cluster cells
        mat <- apply(mat, 1, function(z)
            tapply(z, cut(x, what="species"), mean))
        mat <- apply(mat, 1, function(z)
            tapply(z, cut(x, what="quadrat"), mean))
        colnames(mat) <- spcl$labels
    }
    rownames(mat) <- qcl$labels
    tabasco(mat, use=qcl, sp.ind = spcl, Rowv = FALSE, Colv = FALSE, ...)
    invisible(mat)
}
