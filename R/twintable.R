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
#' @param goodspecies Select \dQuote{good species} for
#'     tabulation. These are either species that were used as
#'     indicator pseudospecies (\code{"indicator"}), or most abundant
#'     species in each final species group breaking ties with
#'     frequency (\code{"leading"}), or \code{"both"} (default). See
#'     \code{\link{goodspec}}. The abundance is estimated after
#'     pseudospecies transformation for all quadrats
#'     and cannot be used together with \code{maxspp}.
#' @param select Select a subset of quadrats.
#'
#' @importFrom vegan vegemite
#'
#' @export
`twintable` <-
    function(object, maxspp, goodspecies, select)
{
    i <- object$quadrat$index
    j <- object$species$index
    ilab <- object$quadrat$labels
    jlab <- object$species$labels
    iclass <- object$quadrat$iclass
    jclass <- object$species$iclass
    mat <- twin2mat(object)
    ## select first a subset
    if (!missing(select)) {
        i <- i[select[i]]
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
