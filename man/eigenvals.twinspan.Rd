% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eigenvals.twinspan.R
\name{eigenvals.twinspan}
\alias{eigenvals.twinspan}
\alias{eigenvals}
\title{Eigenvalues of twinspan Divisions}
\usage{
\method{eigenvals}{twinspan}(x, what = c("quadrat", "species"), ...)
}
\arguments{
\item{x}{\code{\link{twinspan}} result object.}

\item{what}{Return eigenvalues of \code{"quadrat"} or
\code{"species"} divisions.}

\item{\dots}{Other arguments (ignored).}
}
\value{
Vector of eigenvalues ordered by division number, with zero
    for divisions that were not evaluated (terminal groups,
    division terminated in earlier steps).
}
\description{
Function returns the eigenvalues of \code{\link{twinspan}} divisions.
}
\details{
The eigenvalues are for the first correspondence analysis axis of
downweighted pseudospecies data (see \code{\link{twinsform}}). The
eigenvalues are not evaluated for final groups which are not
divided further, nor if the analysis had terminated in a mother
class. This leaves zeros in the vector of eigenvalues.

The eigenvalues are not additive and each is based on slightly
different data due to downweighting (see \code{\link{twinsform}})
and all come from the first axis of respective analysis. They may
have a relation to the magnitude of differences in division, but
there is no information on the total Chi-square nor on its
reduction in division: the eigenvalues only describe the strength
of the axis used as a device of division. Neither can eigenvalues
for quadrats and species compared: data for species clustering are
constructed in a completely different way than data for quadrat
classification.

Eigenvalues are in general not decreasing, but they can increase
when divisions proceed. Function
\code{\link{as.dendrogram.twinspan}} can try to use division
eigenvalues as dendrogram heights, but the resulting trees often
have inversions and look messy or even unreadable.
}
\examples{

data(ahti)
tw <- twinspan(ahti)
eigenvals(tw)

}
\seealso{
\code{\link[vegan]{eigenvals}} in \CRANpkg{vegan}.
    \code{\link{summary.twinspan}} displays the eigenvalues
    numerically, and they can be used in
    \code{\link{as.dendrogram.twinspan}}.
}
