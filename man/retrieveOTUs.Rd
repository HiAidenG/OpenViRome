% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseSQL.R
\name{retrieveOTUs}
\alias{retrieveOTUs}
\title{retrieveOTUs}
\usage{
retrieveOTUs(conn = connectToDatabase(), ...)
}
\arguments{
\item{conn}{A connection to the Serratus SQL database. If not provided, a
connection will be established.}

\item{...}{A list of filters to apply to the query. Filters are provided as
named arguments. For example, to filter by runID, you would provide
runID = "SRR1234567". If a filter's value is NA (i.e. runID = NA),
the filter will be ignored, but the virome will be returned with
the specified column.}
}
\value{
A tibble containing the requested OTUs and columns. Note that more
rows than sOTUs may be returned if an sOTU is associated with
multiple runs in the data requested.
}
\description{
Retrieves OTUs from the Serratus SQL based on the provided
filters. This function is the main entry point for querying the database.
As of the current version, the function supports the following filters:
\itemize{
\item runID: SRA run accessions
\item biosampleID: biosample accessions
\item bioprojectID: bioproject accessions
\item librarySource: metadata-annotated source organism for the library. Note
that you may provide a taxon at any level. openviRome will
automatically collect child taxa.
\item biosampleTissue: metadata-annotated tissue from which the biosample was
derived.
\item phylum: Predicted taxonomic phylum of the sOTU.
\item family: Predicted family of the sOTU.
\item genus: Predicted genus of the sOTU.
\item order: Predicted order of the sOTU.
\item species: Predicted species of the sOTU.
\item normalizedNodeCoverage: Minimum normalized node coverage of the sOTU
(within the run). Normalization is performed by
dividing the coverage for that sOTU by total
spots in the run, scaling by 1e6.
\item STAT: Filter for runs containing reads mapping to a specified taxon.
\item minKmerPerc: The minimum percentage of kmers in the library for the given
STAT value.
}
}
\examples{
retrieveOTUs(runID = NA, biosampleID = "SAMN12345678")
retrieveOTUs(runID = NA, librarySource = c("Homo sapiens", "Apicomplexa"),
normalizedNodeCoverage = 20)
retrieveOTUs(runID = NA, biosampleID = NA, bioprojectID = NA, librarySource =
"Homo sapiens", biosampleTissue = "heart")

}
