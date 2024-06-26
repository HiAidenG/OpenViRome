% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseSQL.R
\name{handleTaxSearch}
\alias{handleTaxSearch}
\title{handleTaxSearch}
\usage{
handleTaxSearch(taxNames = NULL)
}
\arguments{
\item{taxNames}{A vector of taxonomic names.}
}
\value{
A vector of species names corresponding to the input taxa.
}
\description{
Given a set of taxonomic names, collect all species belonging
to those taxa, returning them as names.
}
\examples{
handleTaxSearch(taxNames = c("Nematoda", "Arthropoda")))

}
\keyword{internal}
