% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseSQL.R
\name{getRunDiff}
\alias{getRunDiff}
\title{getRunDiff}
\usage{
getRunDiff(allData = NULL, virome = NULL)
}
\arguments{
\item{allData}{A data frame containing all data processed by Serratus for
a given query. See getAllDataProcessed. MUST contain a 'runID' column.}

\item{virome}{A data frame containing the virome data. MUST contain a
'runID' column.}
}
\value{
A data frame containing the rows in allData that are not present in
virome.
}
\description{
Returns the rows in allData that are not present in virome.
Intended for use with the getAllDataProcessed function.
}
\examples{
conn <- connectToDatabase()
allData <- getAllDataProcessed(conn, runID = NA, librarySource = "Nematoda")
virome <- getViromeData(conn, runID = NA, librarySource = "Nematoda")
getRunDiff(allData, virome)
}
