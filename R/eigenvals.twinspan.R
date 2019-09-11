#' Eigenvalues of twinspan Divisions
#'
#' @param x \code{\link{twinspan}} result object.
#' @param what Return eigenvalues of \code{"quadrat"} or
#'     \code{"species"} divisions.
#' @param \dots Other arguments (ignored).
#' @importFrom vegan eigenvals
#' @export eigenvals
#' @export
`eigenvals.twinspan` <-
    function(x, what = c("quadrat", "species"), ...)
{
    what <- match.arg(what)
    out <- x[[what]]$eig
    names(out) <- paste0("Div", seq_along(out))
    out <- out[out > 0]
    attr(out, "sumev") <- NA
    class(out) <- "eigenvals"
    out
}
