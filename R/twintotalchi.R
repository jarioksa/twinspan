#' Total Scaled Chi-square of All Divisions and Groups
#'
#' Total scaled Chi-square is the same as the sum of all eigenvalues
#' in Correspondence Analysis. It is a measure of heterogeneity of a
#' set of data. The function finds this measure for all divisions and
#' terminal groups in \code{\link{twinspan}}
#'
#' Function first reconstructs the data as it was internally used in
#' \code{\link{twinspan}} using functions \code{\link{twin2stack}} and
#' \code{\link{twin2specstack}}. The scaled Chi-square is calculated
#' with support function \code{totalchi} that can be called
#' independently for any matrix with non-negative data. The scaling in
#' Chi-square means that data are standardized to unit sum before
#' calculating the actual Chi-square. This is often called the sum of
#' all eigenvalues in Correspondence Analysis, but no eigenvalues are
#' evaluated in these functions, only their potential sum.
#'
#' @param x A matrix of non-negative data for \code{totalchi} or a
#'     \code{\link{twinspan}} result object for \code{twintotalchi}.
#' @param what Analyse \code{quadrat} or \code{species}
#'     classification.
#'
#' @seealso The basic functions are \code{\link{twin2stack}} and
#'     \code{\link{twin2specstack}} that construct the data matrices.
#'
#' @examples
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' twintotalchi(tw)

### sum of all eigenvalues in Correspondence Analysis (total Chi-squared)
#' @rdname twintotalchi
#' @export
`totalchi` <-
    function(x)
{
    x <- x/sum(x)
    rc <- outer(rowSums(x), colSums(x))
    x <- (x - rc)/sqrt(rc)
    sum(x^2)
}

### Evaluates the sum of all eigenvalues of Correspondence Analysis
### for all divisions and terminal groups > 1 member. Twinspan only
### evaluated the first eigenvalue of divisions

#' @rdname twintotalchi
#' @export
`twintotalchi` <-
    function(x, what = c("quadrat", "species"))
{
    what <- match.arg(what)
    chi <- numeric(2^(x$levelmax+1) - 1)
    for(lev in 0:x$levelmax) {
        ids <- cut(x, level = lev, what = what)
        tab <- table(ids)
        for (k in unique(ids))
            if(sum(ids==k) > 1) {
                z <- switch(what,
                    "quadrat" =
                        twin2stack(x, subset = ids==k, downweight = TRUE),
                    "species" =
                        twin2specstack(x, subset = ids==k, downweight = TRUE))
                chi[k] <- totalchi(z)
            }
    }
    chi
}

