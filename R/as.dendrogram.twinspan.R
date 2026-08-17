#' Extract Species or Quadrat Dendrograms
#'
#' Function extracts the species or quadrat classification as a
#' hierarchic \code{\link[stats]{dendrogram}}. Unlike \code{hclust},
#' twinspan dendrograms show all quadrats or species, and the final
#' divisions are polytomous.
#'
#' The dendrogram heights are levels of divisions, or total
#' Chi-squares of divisions and groups, or eigenvalues of divisions
#' depending on argument \code{height}. Terminal groups have no
#' eigenvalues, because they were not considered for division. For
#' them the method uses arbitrary value that for a group of \eqn{n}
#' units is proportion \eqn{(n-1)/n} of the height of mother
#' division. Chi-squares are evaluated also for terminal groups.
#' There is no guarantee that eigenvalues or Chi-squares decrease in
#' divisions, and there may be reversals where lower levels are higher
#' than their mother groups, and the plotted trees can be messy and
#' unreadable. Chi-squares decrease more monotonically than
#' eigenvalues of first axis.
#'
#' Function allows selecting a subset of items.
#'
#' \R{} has a wealth of functions to handle and display
#' dendrograms. See \code{\link[stats]{dendrogram}} for general
#' description. There is even stronger support in packages (for
#' instance, \CRANpkg{dendextend}).
#'
#' The terminal groups of \code{\link{twinspan}} trees are not binary,
#' but may have several elements (quadrats, species).  In
#' \code{\link[stats]{dendrogram}} plots, it is best to set
#' \code{type="triangle"} for nicer looking trees.
#'
#' @return A \code{\link[stats]{dendrogram}} object.
#'
#' @examples
#'
#' ## Large datasets are difficult to show in dendrograms, but
#' #'subset' allows slicing dendrogram.
#'
#' data(ahti)
#' tw <- twinspan(ahti)
#' ## cut at level 2 into groups
#' cl2 <- cut(tw, level=2)
#' table(cl2)
#' den <- as.dendrogram(tw, subset = cl2 == 4)
#' str(den, max.level = 4)
#' plot(den, type = "triangle")
#' den <- as.dendrogram(tw, height="chi", subset = cl2 == 4)
#' plot(den, type = "triangle")
#' ## show only most frequent species
#' freq <- colSums(ahti > 0)
#' den <- as.dendrogram(tw, what = "species", height="chi", subset = freq > 16)
#' plot(den, type = "triangle")
#'
#' @seealso \code{\link{as.hclust.twinspan}}, \code{\link{dendrogram}}.
#'
#' @param object \code{\link{twinspan}} result object.
#' @param height Use either division levels (\code{"level"}), total
#'     Chi-squares (\code{"chi"}) or eigenvalues of first axis
#'     (\code{"eigen"}) of division as dendrogram heights.
#' @param subset A logical vector or indices that select a subset of
#'     items to a dendrogram.
#' @param what Return either a \code{"quadrat"} or \code{"species"}
#'     dendrogram.
#'
#' @param \dots Other parameters to functions.
#'
#' @importFrom stats as.dendrogram
#'
#' @export
`as.dendrogram.twinspan` <-
    function(object, height = c("level", "chi", "eigen"),
             what = c("quadrat", "species"), subset, ...)
{
    what <- match.arg(what)
    height <- match.arg(height)
    obj <- object[[what]]
    clid <- cut(object, what=what)
    if (!missing(subset)) {
        clid <- clid[subset]
        obj$labels <- obj$labels[subset]
    }
    len <- length(obj$eig) * 2 + 1
    state <- character(len)
    state[unique(clid)] <- "leaf"
    for(k in rev(seq_len(len %/% 2)))
        if(state[2*k] != "" || state[2*k+1] != "")
            state[k] <- "branch"
    ## eigen: divisions have eigenvalue, for terminal group of size n
    ## use proportion (n-1)/n of mother division
    if (height == "eigen") {
        eig <- numeric(len)
        eig[seq_along(obj$eig)] <- obj$eig
        for(k in which(state == "leaf")) {
            nk <- sum(clid == k)
            eig[k] <- (1 - 1/nk) * eig[k %/% 2]
        }
    } else if (height == "chi") {
        eig <- twintotalchi(object, what = what)
    } else {
        pow2 <- 2^(0:(object$levelmax+1))
        hmax <- sum(max(which(nchar(state) > 0 )) >= pow2) + 1
    }
    z <- list()
    ## In the beginning state can be "", "leaf" or "branch". At the
    ## end of a cycle, processed item is re-labelled "done".
    for(k in rev(seq_along(state))) {
        if(state[k] == "")
            next
        if(state[k] == "leaf") {
            zk <- as.list(which(clid == k))
            attr(zk, "members") <- length(zk)
            attr(zk, "midpoint") <- (length(zk)-1)/2
            labs <- obj$labels[clid == k]
            if (height %in% c("chi","eigen")) {
                hi <- 0
            } else {
                hi <- hmax - sum(k >= pow2) - 1
            }
            for (i in seq_len(length(zk))) {
                attr(zk[[i]], "label") <- labs[i]
                attr(zk[[i]], "members") <- 1L
                attr(zk[[i]], "height") <- hi
                attr(zk[[i]], "leaf") <- TRUE
            }
        }
        else { # a branch
            k1 <- 2*k
            k2 <- 2*k+1
            ## with subset a branch may have lost one of its children
            if (state[k1] == "" || state[k2] == "")
                next
            ## child can be an undone branch: go deeper
            repeat{
                if (state[k1] == "")     # follow other branch
                    k1 <- k1 + 1
                if (state[k1] == "done") # done!
                    break
                k1 <- 2 * k1             # k1 was "branch": go deeper
            }
            repeat{
                if (state[k2] == "")
                    k2 <- k2 + 1
                if (state[k2] == "done")
                    break
                k2 <- 2 * k2
            }

            x <- c(k1, k2)
            x <- as.character(x)
            zk <- list(z[[x[1]]], z[[x[2]]])
            attr(zk, "members") <- attr(z[[x[1]]], "members") +
                attr(z[[x[2]]], "members")
            attr(zk, "midpoint") <- (attr(z[[x[1]]], "members") +
                                     attr(z[[x[1]]], "midpoint") +
                                     attr(z[[x[2]]], "midpoint"))/2
            z[[x[1]]] <- z[[x[2]]] <- NULL
        }
        ## Divisions have eigenvalue, but ev is never evaluated for
        ## terminal groups ("leaf"). We use an arbitrary value: for
        ## group of size n proportion (n-1)/n of the eigenvalue of
        ## mother division.
        if (height %in% c("chi","eigen")) {
            attr(zk, "height") <- eig[k]
        }
        else
            attr(zk, "height") <- hmax - sum(k >= pow2)
        z[[as.character(k)]] <- zk
        ## done!
        state[k] <- "done"
    }
    structure(z[[1]], class="dendrogram")
}
