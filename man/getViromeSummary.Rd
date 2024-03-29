% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viromeQueries.R
\name{getViromeSummary}
\alias{getViromeSummary}
\title{getViromeSummary}
\usage{
getViromeSummary(virome = NULL)
}
\arguments{
\item{virome}{A virome object}
}
\value{
A list of summary statistics (# of unique sOTUs, mean normalized
coverage, median normalized coverage, max normalized coverage, virus
positive runs, total runs processed).
}
\description{
Return a number of summary statistics for a given virome.
}
\examples{
con <- palmid::SerratusConnect()
virome <- getVirome(tax="Canis", con=con)
summary <- getViromeSummary(virome = virome)

}
