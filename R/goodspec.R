### Function to select important species to display. These are either
### (1) species used as indicator species, or (2) most abundant and
### common species of final units of species classification, or (3) both.

#' Most Diagnostic and Characteristic Species
#'
#' Function can return (1) species that are used as indicators in
#' \code{twinspan} classification, or (2) most abundant and common
#' species (known as leading species) of each final \code{twispan}
#' class, or both. The default is to return indices that can be used
#' in subsetting data, but it is also possible to return the species
#' labels.
#'
#' @return Either indices of species or their name labels.
#'
#' @seealso Functions \code{\link{image.twinspan}},
#'     \code{\link{twinheat}} and \code{\link{twintable}} can use the
#'     function to select species.
#'
#' @examples
#' ## similar to plot(tw, what="species") but labels classes by their
#' ## leading species
#' data(ahti)
#' tw <- twinspan(ahti)
#' plot(as.dendrogram(tw, what = "species", subset = goodspec(tw, "leading")))
#'
#' @param object \code{\link{twinspan}} result object.
#' @param what Select either species that are used as
#'     \code{"indicator"} pseudospecies or most abundant and common
#'     (or \code{"leading"}) species of each group species
#'     classification or \code{"both"}.
#' @param label Return species labels instead of indices.
#'
#' @export
`goodspec` <-
    function(object, what = c("both", "indicator", "leading"), label = FALSE)
{
    what <- match.arg(what)
    ## pick indicator species
    if (what %in% c("both", "indicator")) {
        inds <- abs(object$quadrat$indicators)
        inds <- inds[inds > 0]
        inds <- object$quadrat$pseudo2species[inds]
        inds <- sort(unique(inds))
    }
    ## pick most abundant species (and try breaking ties with
    ## frequency) for each final group of species classification
    if (what %in% c("both", "leading")) {
        mat <- twin2mat(object)
        abu <- colSums(mat)
        cnt <- colSums(mat > 0)
        cl <- cut(object, what = "species")
        id <- sort(unique(cl))
        lead <- numeric(length(id))
        for (k in seq_along(id)) {
            i <- which(cl==id[k])
            lead[k] <- i[order(abu[i], cnt[i], decreasing=TRUE)[1]]
        }
        ## combine with indicators if calculated
        if (what == "both") {
            inds <- sort(unique(c(lead,inds)))
        } else { # sorted by class number
            inds <- lead
        }
    }
    if (label)
        inds <- object$species$labels[inds]
    inds
}

