# data.R

#' TylenchoideaVirome
#'
#' Filtered Serratus SQL data with the member genus of the Tylenchoidea
#' superfamily.
#'
#' @format A list of data frames. The first has 611 rows and 20 columns.
#' Each row represents a single palm_id, with the corresponding run in which it was identified, as well as the
#' viral species to which is belongs (sOTU). Columns as follows:
#' \describe{
#'  \item{run}{SRA run accession number.}
#'  \item{scientific_name}{Metadata annotation for the library source species.}
#'  \item{bio_sample}{SRA BioSample accession number.}
#'  \item{bio_project}{Source BioProject accession number.}
#'  \item{palm_id}{Serratus identifier for this palm print.}
#'  \item{sotu}{Serratus identifier for the clustered sOTU (species-like
#'  operational taxonomic unit).}
#'  \item{nickname}{Serratus nickname for the sOTU.}
#'  \item{gb_acc}{GenBank accession for the closest aligned sequence via BLAST}
#'  \item{gb_pid}{Percent identity of the BLAST alignment.}
#'  \item{gb_eval}{E-value for the BLAST alignment.}
#'  \item{tax_species}{Taxonomic species annotation for the GenBank aligned
#'  virus.}
#'  \item{node}{Node number for the vOTU.}
#'  \item{node_coverage}{Coverage of the vOTU in this library.}
#'  \item{node_pid}{Node percent identity.}
#'  \item{node_eval}{Node E-value.}
#'  \item{node_qc}{Bool; whether the node passed quality control checks.}
#'  \item{node_seq}{Detected palm sequence for this vOTU.}
#'  \item{node_coverage_norm}{Normalized coverage of the vOTU in this library.}
#'  \item{tax_phylum}{Taxonomic phylum annotation for the GenBank aligned sequence}
#'  }
#' The second data frame has 408 rows and 5 columns. This is a table of all
#' runs queried by Serratus annotated as being associated with the 'Tylenchoidea'
#' superfamily. Columns as follows:
#' \describe{
#' \item{run}{SRA run accession number.}
#' \item{scientific_name}{Metadata annotation for the library source species.}
#' \item{tax_id}{NCBI taxonomy identifier for the library source species.}
#' \item{spots}{Number of spots (reads) in the library.}
#' \item{virus_positive}{Bool; whether a palm print was detected in this library.}
#'
#' @references
#' Edgar, R.C., Taylor, B., Lin, V. et al. Petabase-scale sequence alignment
#' catalyses viral discovery. Nature 602, 142–147 (2022).
#' https://doi.org/10.1038/s41586-021-04332-2
#' @examples
#' \dontrun {
#' TylenchoideaVirome
#' }
#'
"TylenchoideaVirome"

# [END]
