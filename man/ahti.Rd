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
to use \code{\link{browseVignettes}} (but see Examples below).

Oksanen & Ahti (1982) allocated the quadrats into vegetation classes
that they called \sQuote{noda}. The name of quadrat is based on the
original nodum followed by the original number of the quadrat. In
northern Finland (nortern boreal zone, roughly corresponding to the are
of intensive reindeer grazing) five noda were separated: \code{Ster} for
\emph{Stereocaulon}, \code{Vuli} for \emph{Vaccinium uliginosum},
\code{Vmyr} for \emph{Vaccinium myrtillus}, \code{Call} for
\emph{Calluna}, and \code{Empe} for \emph{Empetrum}. No noda were
separated in middle boreal (names starting with \code{M}) and southern
boreal (names starting with \code{S}) zones, but these were divided to
young (\code{Myng}, \code{Syng}) and old (\code{Mold}, \code{Sold})
forests. The only exception was that the vegetation on southerly slopes
of eskers was allocated to \emph{Thymus} nodum (\code{Thym}). Hemiboreal
forests on the southwest coast were allocated to nodum \code{Hbor} with
no age distinction.

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
