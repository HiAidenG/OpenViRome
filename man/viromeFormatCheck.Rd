% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viromeQueries.R
\name{viromeFormatCheck}
\alias{viromeFormatCheck}
\title{viromeFormatCheck}
\usage{
viromeFormatCheck(virome = NULL)
}
\arguments{
\item{virome}{A virome object}
}
\value{
A boolean value
}
\description{
Checks if the virome object is in the correct format. Mostly for
use with the shiny frontend.
}
\examples{
con <- palmid::SerratusConnect()
virome <- getVirome(tax="Canis", con=con)
if (viromeFormatCheck(virome = virome)) {
 print("Virome is in the correct format.")
}

}
