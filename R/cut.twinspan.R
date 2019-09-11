#' Return twinspan Classification at Given Level
#'
#' Function returns a vector of \code{twinspan} classes at given level
#' of hierarcy for quadrats.
#'
#' @param x \code{twinspan} result.
#' @param level Level of hierarchy for classification. If missing, the
#'     highest lvel used in the object will be returned.
#' @param what Return either a \code{"quadrat"} or \code{"species"}
#'     classification vector.
#' @param \dots Other parameters (ignored).
#'
#' @return A vector of class numbers for the given level of hierarchy
#'     using Twinspan identifiers. A mother class can be found through
#'     integer division by 2 for a given class identifier.
#' @export
`cut.twinspan` <-
    function(x, level, what = c("quadrat", "species"), ...)
{
    what <- match.arg(what)
    if (missing(level))
        level <- 31L
    clmax <- 2^(level+1) - 1
    cl <- x[[what]]$iclass
    big <- cl > clmax
    while(any(big)) {
        ## mother class by integer division
        cl[big] <- cl[big] %/% 2
        big <- cl > clmax
    }
    cl
}
