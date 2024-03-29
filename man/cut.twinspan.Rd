% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cut.twinspan.R
\encoding{UTF-8}
\name{cut.twinspan}
\alias{cut.twinspan}
\alias{cuth}
\title{Return twinspan Classification at Given Level}
\usage{
\method{cut}{twinspan}(x, level, what = c("quadrat", "species"), binname = FALSE, ...)

cuth(x, what = c("quadrat", "species"), ngroups, binname = FALSE)
}
\arguments{
\item{x}{\code{twinspan} result.}

\item{level}{Level of hierarchy for classification. If missing, the
final level used in the object will be returned.}

\item{what}{Return either a \code{"quadrat"} or \code{"species"}
classification vector.}

\item{binname}{Use binary label for classes instead of decimal number.}

\item{\dots}{Other parameters (ignored).}

\item{ngroups}{Number of groups.}
}
\value{
A vector of class numbers for the given level of hierarchy
    using Twinspan identifiers. For identifiers and levels, see
    \code{\link{summary.twinspan}},
    \code{\link{as.hclust.twinspan}}.
}
\description{
Returns a vector of \code{twinspan} classes at a given level of
hierarchy or classes respecting group heterogeneity for quadrats or
species.
}
\details{
\code{\link{twinspan}} returns only the the classification at the
final level, but any upper level classes can be found by integer
divisions by 2. \code{\link{twinspan}} bases classification
principally on splitting polished ordination axis. Sometimes this
allocation is in conflict with indicator pseudospecies. This is
called \sQuote{misclassification} (see
\code{\link{misclassified}}). Function
\code{\link{predict.twinspan}} returns the similar classification
based on indicator pseudospecies (also for new data), and
\code{\link{misclassified}} analyses the differences of these
classifications.

Function \code{cuth} cuts the classification by class heterogeneity
instead of level, and can be used to implement the modified method
of Roleček et al. (2009). The groups are formed with decreasing
heterogeneity but respecting the hierarchy. Total chi-square (also
known as inertia) is used as the criterion of heterogeneity. The
criterion is calculated with \code{\link{twintotalchi}} and the
criterion is based on the same data matrix as internally used in
\code{\link{twinspan}}. The function can also be used for species
classification, also with the internally used modified species
matrix.
}
\examples{

data(ahti)
tw <- twinspan(ahti)
cut(tw)
## traditional twinspan classification by level of hierarchy
cut(tw, level=3)
cut(tw, what = "species")
## number of groups as with level=3, but by group heterogeneity
cuth(tw, ngroups = 8)

}
\references{
Roleček, J, Tichý, L., Zelený, D. & Chytrý, M. (2009). Modified
TWINSPAN classification in which the hierarchy respects cluster
heterogeneity. \emph{J Veg Sci} 20: 596--602.
}
\seealso{
\code{\link{predict.twinspan}} gives similar classes, but
    based on indicator pseudospecies. \code{\link{cutree}} provides
    a similar functionality for \code{\link{hclust}}
    trees. Function \code{\link{as.hclust.twinspan}} generates
    corresponding tree presentation, and
    \code{\link{plot.twinspan}} will print that tree labelling
    internal nodes (divisions).
}
