\name{ahti}
\alias{ahti}
\docType{data}
\title{
Vegetation of Lichen-rich Pine Forests in Finland
}

\description{
Vegetation data set of 170 quadrats and 115 species with visually
estimated cover percentage values (Oksanen & Ahti, 1982).
}

\usage{data("ahti")}

\details{

The taxon names are given with 4+4 acronyms of the binomial name. The
complete scientific names are listed in \code{\link{vignette}}
\code{"ahtinames"}; perhaps the easiest way of reading the vignette is
to use \code{\link{browseVignettes}} (but see Examples below). Only
understorey plants are included in the current file: trees and shrubs,
also as seedlings, are not included, although they are listed in the
original article.

Oksanen & Ahti (1982) allocated the quadrats into vegetation classes
that they called \sQuote{noda}. The name of quadrat is based on the
original nodum followed by the original number of the quadrat. In
northern Finland (nortern boreal zone, roughly corresponding to the
are of intensive reindeer grazing) five noda were separated:
\code{Ster} for \emph{Stereocaulon}, \code{Vuli} for \emph{Vaccinium
uliginosum}, \code{Vmyr} for \emph{Vaccinium myrtillus}, \code{Call}
for \emph{Calluna}, and \code{Empe} for \emph{Empetrum}. Most Middle
boreal, south boreal and hemiboreal stands were allocated into
\emph{Cetraria islandica} nodum, but tabulated separately for young
and old forest stands and by vegetation zone: middle boreal (names
starting with \code{M}) and southern boreal (names starting with
\code{S}) zones divided into young (\code{Myng}, \code{Syng}) and old
(\code{Mold}, \code{Sold}) forests. Hemiboreal forests (\code{Hbor})
were not divided by age classes. The only separate southern boreal
type was \emph{Thymus serpyllum} nodum (\code{Thym}) on the sunny
slopes of eskers.

}

\source{
Oksanen & Ahti (1982).
}

\references{
Oksanen, J. & Ahti, T. (1982) Lichen-rich pine forest vegetation in
Finland. \emph{Annales Botanici Fennici} \bold{19,} 275--301.
}

\examples{
## Read the vignette listing complete scientific names
if (interactive()) {
vignette("ahtinames", package="twinspan")
}
}

\keyword{datasets}
