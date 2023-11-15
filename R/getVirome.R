#' @title getVirome
#' @description Returns a virome object given either:
#' 1) a taxonomic term and (optional) taxonomic level
#' 2) a character vector of SRA accessions
#' Note: the broader the taxonomic search term, the longer the query will
#' take.
#' @param tax A taxonomic term
#' @param sra A character vector of SRA accessions
#' @return A virome object
#' @export
#' @examples
#' con <- palmid::SerratusConnect()
#' getVirome(tax = "Meloidogyne", con = con)
#' getVirome(sra = c("SRR17756040", "SRR5942326"), con = con)
#' @importFrom dplyr tbl collect filter mutate group_by summarise n
#' @importFrom taxize classification downstream
#' @references
#' Wickham H, François R, Henry L, Müller K, Vaughan D (2023).
#' dplyr: A Grammar of Data Manipulation.
#' https://dplyr.tidyverse.org, https://github.com/tidyverse/dplyr.
#'
#' A. Babaian and R. C. Edgar (2022),
#'Ribovirus classification by a polymerase barcode sequence. PeerJ
#'
#' R. C. Edgar et al. (2021),
#' Petabase-scale sequence alignment catalyses viral discovery. Nature.
#' @export
getVirome <- function(tax = NULL, sra = NULL, con = NULL) {

  if (is.null(con)) {
    stop("Please provide a connection to the Serratus database
         (see palmid::SerratusConnect)")
  }
  if (is.null(tax) & is.null(sra)) {
    stop("Must provide either a taxonomic term or a character vector of
         SRA accessions")
  }
  if (!is.null(tax) & !is.null(sra)) {
    stop("Must provide either a taxonomic term or a character vector of
         SRA accessions")
  }
  if (!is.null(tax) & is.null(sra)) {

    # Figure out taxonomic level
    tax <- stringr::str_to_title(tax)

    # Check if a genus was provided
    rank <- taxize::classification(tax, db = 'ncbi')[[1]]
    rank <- rank[rank$name == tax, 'rank']
    if (is.null(rank)) {
      stop("Error: could not find taxonomic term in NCBI taxonomy database")
    }
    else if (rank == 'genus') {
      searchTerms <- tax
    }
    else {
      # Have to go down to genus level
      searchTerms <- taxize::downstream(tax, db ='ncbi', downto = 'genus')
      searchTerms <- searchTerms[[1]]['childtaxa_name'][[1]]
    }

    # Get the virome object
    regexPattern <- paste(searchTerms, collapse = "|")
    virome <- tbl(con, "palm_virome") %>%
      dplyr::filter(grepl(pattern = regexPattern, x = scientific_name,
                          ignore.case = TRUE)) %>%
      dplyr::collect()

    return(virome)
  }

  else if (is.null(tax) & !is.null(sra)) {
    # Get the virome object
    virome <- tbl(con, "palm_virome") %>%
      dplyr::filter(run %in% sra) %>%
      dplyr::collect()

    return(virome)

  }
}
