### TWINSPAN table

#' Community Table Ordered by Twinspan Classification
#'
#' Prints a community table of pseudospecies ordered by
#' \code{\link{twinspan}} classification.
#'
#' Function prints a compact community table of pseudospecies
#' values. The table is ordered by clustering both species and
#' quadrats similarly as in \code{\link{summary.twinspan}} or in plot
#' of \code{\link{as.dendrogram.twinspan}}. The classification of each
#' quadrat and species is shown by a sequence of \code{0} and \code{1}
#' indicating division of each level. This string is binary
#' presentation of the decimal class number without the leading
#' \code{1}.
#'
#' Only one character is used for each abundance, and the table is
#' very compact. However, large tables can be divided over several
#' pages or screen windows. The width of the displayed table is
#' controlled by \R{} option \code{width} (see
#' \code{\link{options}}). It is possible to select only a
#' \code{subset} of the quadrats for tabulation giving narrower
#' tables. The number of species can be reduced by setting the maximum
#' number of most abundant species, or alternatively, by restricting
#' tabulation only to \dQuote{good species} which are the most
#' abundant species of each species group (ties broken by species
#' frequency), or species used as indicators, or both.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' ## complete table would be large, but we subset
#' twintable(tw, subset = cut(tw, 2) == 4, goodspecies = "both")
#'
#' @param object \code{\link{twinspan}} result object.
#' @param maxspp Maximum number of most abundant species
#'     displayed. The abundance is estimated with pseudospecies cut
#'     levels. The default is to show all species.
#' @param goodspecies Select \dQuote{good species} for
#'     tabulation. These are either species that were used as
#'     indicator pseudospecies (\code{"indicator"}), or most abundant
#'     species in each final species group breaking ties with
#'     frequency (\code{"leading"}), or \code{"both"} (default). The
#'     abundance is estimated after pseudospecies transformation for
#'     all quadrats and cannot be used together with \code{maxspp}.
#' @param subset Select a subset of quadrats.
#'
#' @importFrom vegan vegemite
#'
#' @export
`twintable` <-
    function(object, maxspp, goodspecies, subset)
{
    i <- object$quadrat$index
    j <- object$species$index
    ilab <- object$quadrat$labels
    jlab <- object$species$labels
    iclass <- object$quadrat$iclass
    jclass <- object$species$iclass
    mat <- twin2mat(object)
    ## select first a subset
    if (!missing(subset)) {
        i <- i[subset[i]]
    }
    ## First see if we want to have only the "good species"
    if (!missing(goodspecies)) {
        if (!missing(maxspp))
            stop("maxspp and goospecies cannot be defined together")
        gd <- goodspec(object, goodspecies)
        j <- j[j %in% gd]
    }
    ## take only maxspp most abundant species - but if there are ties,
    ## take all tied species, even if we go over maxspp
    if (!missing(maxspp) && ncol(mat) > maxspp) {
        sptot <- colSums(mat[i,,drop=FALSE])
        maxlim <- sort(sptot, decreasing = TRUE)[maxspp]
        j <- j[sptot[j] >= maxlim]
    }
    ## add classification strings to names
    jnam <- addbin2name(jclass[j], jlab[j])
    inam <- addbin2name(iclass[i], ilab[i])
    mat <- mat[i,j]
    dimnames(mat) <- list(inam, jnam)
    vegemite(mat, zero="-")
}

### Unexported function to turn the class number into binary string.
### This is the same as the class number as binary number, but strips
### the leading '1'.
`class2bin` <-
    function(cl)
{
    str <- ""
    while(cl > 1) {  # cl > 0 would give the binary number
        str <- paste0(cl %% 2L, str)
        cl <- cl %/% 2L
    }
    str
}

### Unexported function that adds the binary classification string to
### the item name making all strings equally long

`addbin2name` <-
    function(class, name)
{
    bin <- character(length(class))
    for(i in seq_along(bin)) {
        bin[i] <- class2bin(class[i])
    }
    strlens <- nchar(bin) + nchar(name)
    mxlen <- max(strlens) + 1
    pad <- character(length(class))
    for (i in seq_along(pad))
        pad[i] <- paste0(rep(" ", mxlen-strlens[i]), collapse="")
    paste0(bin, pad, name)
}
