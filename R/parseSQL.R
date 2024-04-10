# Path: OpenViRome/R/parseSQL.R


#' @title retrieveOTUs
#' 
#' @description Retrieves OTUs from the Serratus SQL based on the provided 
#' filters. This function is the main entry point for querying the database.
#' As of the current version, the function supports the following filters:
#'   - runID: SRA run accessions
#'   - biosampleID: biosample accessions
#'   - bioprojectID: bioproject accessions
#'   - librarySource: metadata-annotated source organism for the library. Note
#'                    that you may provide a taxon at any level. openviRome will
#'                    automatically collect child taxa.
#'   - biosampleTissue: metadata-annotated tissue from which the biosample was 
#'                      derived.
#'   - phylum: Predicted taxonomic phylum of the sOTU.
#'   - family: Predicted family of the sOTU.
#'   - genus: Predicted genus of the sOTU.
#'   - order: Predicted order of the sOTU.
#'   - species: Predicted species of the sOTU.
#'   - normalizedNodeCoverage: Minimum normalized node coverage of the sOTU 
#'                             (within the run). Normalization is performed by
#'                             dividing the coverage for that sOTU by total 
#'                             spots in the run, scaling by 1e6.
#'   - STAT: Filter for runs containing reads mapping to a specified taxon.
#'   - minKmerPerc: The minimum percentage of kmers in the library for the given
#'                  STAT value.
#
#' @param conn A connection to the Serratus SQL database. If not provided, a
#'            connection will be established. 
#' @param ... A list of filters to apply to the query. Filters are provided as
#'           named arguments. For example, to filter by runID, you would provide
#'           runID = "SRR1234567". If a filter's value is NA (i.e. runID = NA),
#'           the filter will be ignored, but the virome will be returned with
#'           the specified column. 
#'
#' @return A tibble containing the requested OTUs and columns. Note that more
#'         rows than sOTUs may be returned if an sOTU is associated with 
#'         multiple runs in the data requested.
#' 
#' @examples
#' retrieveOTUs(runID = NA, biosampleID = "SAMN12345678")
#' retrieveOTUs(runID = NA, librarySource = c("Homo sapiens", "Apicomplexa"),
#' normalizedNodeCoverage = 20)
#' retrieveOTUs(runID = NA, biosampleID = NA, bioprojectID = NA, librarySource =
#' "Homo sapiens", biosampleTissue = "heart")
#' 
#' @export
retrieveOTUs <- function(conn = connectToDatabase(), ...) {
  processedArgs <- processUserArgs((list(...)))
  viromeArgs <- processedArgs[[1]]
  specialArgs <- processedArgs[[2]]
  
  # viromeArgs must always have sOTU
  viromeArgs$sotu <- NA

  # Some special parameters REQUIRE certain columns to be present. Add those if
  # not present:
  requiredCols <- getNecessaryCols(specialArgs = specialArgs, currentArgs = viromeArgs)
  viromeArgs <- c(viromeArgs, requiredCols)

  # Special handling for scientific_name column (librarySource) to facilitate
  # searching at any taxonomic level.
  if ('scientific_name' %in% names(viromeArgs)) {
    viromeArgs$scientific_name <- handleTaxSearch(viromeArgs$scientific_name)
  }

  viromeResult <- getVirome(conn = conn, filters = viromeArgs)

  # Apply the special filters
  viromeResult <- specialFilter(conn = conn, virome = viromeResult, filters = specialArgs)

  # All rows should be distinct
  viromeResult <- tibble(distinct(viromeResult))

  # Translate the columns back to the original names
  colnames(viromeResult) <- SQLnames2Rnames(colnames(viromeResult))

  return(viromeResult)
}

getVirome <- function(conn = connectToDatabase(), filters = NULL) {
  table <- "palm_virome"

  query <- constructQuery(filters = filters, SQLtable = table)

  return(executeQuery(conn = conn, query = query))
}

joinAndCleanup <- function(virome = NULL, data = NULL, by = NULL, joining = NULL) {
  # Print missing rows in data
  if (nrow(data) == 0) {
    warning(paste("No data found for", joining))

    # return an empty dataframe
    return(tibble())
  }

  else {
    virome <- dplyr::left_join(virome, data, by = by)

    # Remove any rows with NA values
    virome <- virome[complete.cases(virome), ]

    return(virome)
  }
}

handleNodeCoverageFilter <- function(conn, virome, filterValue) {
  coverageFilters <- list(run = unique(virome$run), spots = NA)
  coverageData <- executeQuery(conn, constructQuery(filters = coverageFilters, SQLtable = "srarun"))

  # If we have the same sotu appear multiple times in one run, we need to sum the node_coverage
  temp <- virome %>%
    group_by(run, sotu) %>%
    summarise(node_coverage = sum(node_coverage), .groups = 'drop')

  temp <- dplyr::left_join(virome, coverageData, by = "run")

  # Calculate the normalized node coverage
  temp$normalizedNodeCoverage <- temp$node_coverage / temp$spots * 1e6

  # Remove the spots and node_coverage columns
  temp <- temp[, !names(temp) %in% c("spots", "node_coverage")]

  if (!is.na(filterValue)) {
    temp <- temp[temp$normalizedNodeCoverage > filterValue, ]
  }

  return(temp)
}

handleBiosampleTissueFilter <- function(conn, virome, filterValue) {
  tissueFilters <- list(biosample_id = unique(virome$bio_sample), tissue = filterValue)

  tissueData <- executeQuery(conn, constructQuery(filters = tissueFilters, SQLtable = "biosample_tissue"))

  # Rename the columns for simple joining
  names(tissueData) <- c("bio_sample", "tissue")

  virome <- joinAndCleanup(virome, tissueData, by = "bio_sample", joining = "biosample tissue")

  return(virome)
  }

handleTaxonomicFilters <- function(conn = NULL, virome = NULL, taxFilters = NULL) {
  taxFilters$sotu <- virome$sotu

  # Construct a query
  taxQuery <- constructQuery(filters = taxFilters, SQLtable = "palm_tax")

  # Execute the query
  taxData <- executeQuery(conn = conn, query = taxQuery)

  virome <- joinAndCleanup(virome, taxData, by = "sotu", joining = "taxonomic data")

  return(virome)
}

handleSTATFilter <- function(conn = NULL, virome = NULL, STATFilters = NULL) {
  if (!is.na(STATFilters$STAT)) {
    # Need to convert the provided STAT taxa to taxids
    STATFilters$STAT <- sapply(STATFilters$STAT, function(x) retrieveTaxIDs(taxNames = x))

    # Go up (or down) to the order level since this is what's currently stored in the table
    STATFilters$STAT <- taxHierarchyTranslate(taxIDs = STATFilters$STAT, level = "order")
  }
  # Construct the query
  if ("minKmerPerc" %in% names(STATFilters)) {
    STATQuery <- constructQuery(SQLtable = "sra_stat",
                                filters = list(run = virome$run,
                                               taxid = STATFilters$STAT,
                                               kmer_perc = STATFilters$minKmerPerc,
                                               name = NA))
  }
  else {
    STATQuery <- constructQuery(SQLtable = "sra_stat", filters = list(run = virome$run, taxid = STATFilters$STAT))
  }

  # Execute the query
  STATData <- executeQuery(conn = conn, query = STATQuery)
  print(STATData)
  names(STATData)[names(STATData) == "name"] <- "STAT_order"
  virome <- dplyr::left_join(virome, STATData, by = "run", relationship = "many-to-many")
  virome <- virome[, !names(virome) %in% "taxid"]

  return(virome)
}

#' @title specialFilter
#' 
#' @description A function to handle special filters that require additional
#' data from the SQL not present in the palm_virome table.
#' 
#' @param conn A connection to the Serratus SQL database.
#' @param virome A dataframe containing the current virome data.
#' @param filters A list of filters to apply to the virome data.
#' 
#' @return A dataframe containing the virome data with the special filters 
#' applied.
#' 
#' @keywords internal
specialFilter <- function(conn = NULL, virome = NULL, filters = NULL) {
  if ('minKmerPerc' %in% names(filters) && !'STAT' %in% names(filters)) {
    # Add an empty STAT filter
    filters$STAT <- NA
  }

  if ('normalizedNodeCoverage' %in% names(filters)) {
    virome <- handleNodeCoverageFilter(conn, virome, filters$normalizedNodeCoverage)
  }

  if ('biosampleTissue' %in% names(filters)) {
    virome <- handleBiosampleTissueFilter(conn, virome, filters$biosampleTissue)
  }

  taxColumns <- c('tax_phylum', 'tax_order', 'tax_genus')
  if (any(taxColumns %in% names(filters))) {
    taxFilters <- filters[names(filters) %in% taxColumns]
    virome <- handleTaxonomicFilters(conn, virome, taxFilters)
  }

  if ('STAT' %in% names(filters)) {
    STATFilters <- if ('kmer_perc' %in% names(filters)) {
      list(STAT = filters$STAT, kmer_perc = filters$kmer_perc)
    } else {
      list(STAT = filters$STAT)
    }
    virome <- handleSTATFilter(conn, virome, STATFilters)
  }
  return(virome)
}

#' @title handleTaxSearch
#' 
#' @description Given a set of taxonomic names, collect all species belonging
#' to those taxa, returning them as names.
#' 
#' @param taxNames A vector of taxonomic names.
#' 
#' @return A vector of species names corresponding to the input taxa.
#' 
#' @examples
#' handleTaxSearch(taxNames = c("Nematoda", "Arthropoda")))
#' 
#' @importFrom taxizedb taxid2name
#' 
#' @keywords internal
handleTaxSearch <- function(taxNames = NULL) {
  taxIDs <- retrieveTaxIDs(taxNames)

  # Translate the taxIDs to the species level
  taxIDs <- taxHierarchyTranslate(taxIDs = taxIDs, level = "species")

  # Get the tax names
  taxNames <- taxizedb::taxid2name(taxIDs, db = "ncbi")

  # Ensure none of the taxNames themselves contain a single quote (i.e. D'Amelio)
  taxNames <- gsub("'", "''", taxNames)

  return(taxNames)
}

#' @title taxHierarchyTranslate
#' 
#' @description Given a set of taxonomic IDs, translate them to a specified 
#' level in the taxonomic hierarchy. Collects all child taxa.
#' 
#' @param taxIDs A vector of taxonomic IDs (must be NCBI taxIDs).
#' 
#' @param level The level in the taxonomic hierarchy to translate the IDs to.
#' 
#' @return A vector of taxonomic IDs at the specified level in the hierarchy.
#' 
#' @examples
#' taxHierarchyTranslate(taxIDs = c(9606, 10090), level = "order")
#' 
#' @importFrom taxizedb taxa_at downstream
#'
#' @keywords internal
taxHierarchyTranslate <- function(taxIDs = NULL, level = NULL) {
  if (is.null(taxIDs) || is.null(level)) {
    stop("Please provide taxIDs and a level.")
  }

  tryCatch({
    IDs <- taxizedb::taxa_at(taxIDs, db = "ncbi", rank = level, warn = TRUE, missing = "lower")
  },
  error=function(cond) {
    stop("Unable to retrieve taxonomic IDs.")
  })

  # Get just the IDs, save it to a vector
  IDs <- unlist(lapply(IDs, function(x) x$id))

  if (length(IDs) == 0) {
    # Try going downstream
    IDs <- taxizedb::downstream(taxIDs, db = "ncbi", downto = level)

    # Get just the IDs, save it to a vector
    IDs <- unlist(lapply(IDs, function(x) x$childtaxa_id))
    names(IDs) <- NULL

    if (length(IDs) == 0) {
      stop("Unable to retrieve taxonomic IDs.")
    }
  }

  return(IDs)
}

#' @title retrieveTaxIDs
#' 
#' @description Retrieves the NCBI taxonomic IDs for a given set of taxonomic 
#' names.
#' 
#' @param taxNames A character vector of taxonomic names.
#' 
#' @return A named list of taxonomic IDs, if successful. Otherwise throws an 
#' error.
#' 
#' @examples
#' taxNames <- c("Homo sapiens")
#' retrieveTaxIDs(taxNames)
#' 
#' @importFrom taxizedb name2taxid
#' 
#' @keywords internal
retrieveTaxIDs <- function(taxNames = NULL) {
  if (is.null(taxNames)) {
    stop("Please provide taxonomic names.")
  }

  tryCatch({
    taxIDs <- taxizedb::name2taxid(taxNames, db = "ncbi", out_type = 'uid')
    names(taxIDs) <- NULL
  },
  error=function(cond) {
    stop("Unable to retrieve taxonomic IDs.")
  })

  return(taxIDs)
}

#' @title processUserArgs
#' 
#' @description Processes user-provided arguments, ensuring they are valid and
#' categorizing them into columns (for database selection) and special arguments
#' (requiring further processing).
#' 
#' @param args A named list of arguments to be processed. The names should
#' correspond to the supported parameter names, and the values are the arguments'
#' values.
#' 
#' @return A list containing two elements: `columns`, a named list of arguments
#' mapped to their SQL column names for database queries, and `specialArgs`,
#' a list of arguments that require special handling.
#' 
#' @examples
#' processUserArgs(list(librarySource = 'Homo sapeins', 
#' STAT = c('Toxoplasma gondii'), 
#' minKmerPerc = 0.5))
#' 
#' @keywords internal
processUserArgs <- function(args = NULL) {
  if (is.null(args) || length(args) == 0) {
    stop("Please provide at least one argument!")
  }
  
  if (!is.list(args) || !all(nzchar(names(args)))) {
    stop("All arguments must be named and non-empty.")
  }
  
  # Retrieve parameters using retrieveParamInfo
  normalParams <- retrieveParamInfo(paramType = 'querying')
  specialParams <- retrieveParamInfo(paramType = 'specialQuerying')
  
  allParamNames <- c(normalParams, specialParams)
  invalidArgs <- setdiff(names(args), allParamNames)
  if (length(invalidArgs) > 0) {
    stop("Invalid argument(s): ", paste(invalidArgs, collapse=", "),
         ". Valid arguments are: ", paste(allParamNames, collapse=", "))
  }
  
  # Separate the normal from special arguments
  columns <- args[names(args) %in% normalParams]
  specialArgs <- args[names(args) %in% specialParams]
  
  # Translate column names to the SQL schema using sqlNames2Rparams
  columns <- setNames(columns, Rnames2SQLnames(names(columns)))
  specialArgs <- setNames(specialArgs, Rnames2SQLnames(names(specialArgs)))
  
  return(list(columns = columns, specialArgs = specialArgs))
}

#' @title Rnames2SQLnames
#' 
#' @description Translates R parameter names to SQL column names.
#' 
#' @param Rnames A character vector of R parameter names.
#' 
#' @return A character vector of SQL column names.
#' 
#' @examples
#' Rnames2SQLnames(c("runID", "biosampleID", "bioprojectID"))
#' 
#' @keywords internal
Rnames2SQLnames <- function(Rnames = NULL) {
  if (is.null(Rnames)) {
    stop("Please provide R names.")
  }
  
  # Get the dictionary of all parameters and their SQL names. Names are
  # the R names. Values are the sql names.
  dict <- getSQLLookupDict()
  
  # Get the SQL names
  sqlNames <- dict[Rnames]
  
  return(sqlNames)
}

#' @title SQLNames2Rnames
#' 
#' @description Translates SQL column names to R parameter names.
#' 
#' @param sqlNames A character vector of SQL column names.
#' 
#' @return A character vector of R parameter names.
#' 
#' @examples
#' SQLNames2Rnames(c("run", "biosample_id", "bioproject_id"))
#' 
#' @keywords internal
SQLnames2Rnames <- function(sqlNames = NULL) {
  if (is.null(sqlNames)) {
    stop("Please provide SQL names.")
  }
  
  dict <- getSQLLookupDict()
  
  # Map each SQL name to its R name directly to maintain order
  Rnames <- sapply(sqlNames, function(name) {
    if (name %in% dict) {
      return(names(dict)[dict == name])
    } else {
      stop(paste("SQL name", name, "not found in dictionary."))
    }
  })
  
  return(Rnames)
}

getCategoricalPlottingParams <- function() {
  c("runID", "biosampleID", "bioprojectID", 
    "librarySource", "biosampleTissue", 
    "phylum", "family", "genus", "order",
    "species", "GBAcc", "STAT")
}

getNumericalPlottingParams <- function() {
  c("normalizedNodeCoverage", "GBpID")
}

getQueryingParams <- function() {
  c("runID", "biosampleID", "bioprojectID", 
    "family",
    "species", "GBAcc",
    "librarySource", "GBpID")
}

getSpecialQueryingParams <- function() {
  c("normalizedNodeCoverage", "biosampleTissue", 
    "phylum", "genus", "order", "STAT",
    "kmerPerc")
}

getInequalityQueryingParams <- function() {
  c("minKmerPerc")
}

getSQLLookupDict <- function() {
  c("runID" = "run", "biosampleID" = "bio_sample", "bioprojectID" = "bio_project",
    "librarySource" = "scientific_name", "biosampleTissue" = "tissue",
    "phylum" = "tax_phylum", "family" = "tax_family", "genus" = "tax_genus", "order" = "tax_order",
    "species" = "tax_species", "GBAcc" = "gb_acc", "STAT" = "STAT",
    "normalizedNodeCoverage" = "normalizedNodeCoverage", "GBpID" = "gb_pid",
    "minKmerPerc" = "kmer_perc", "sOTU" = "sotu")
}

getNecessaryCols <- function(specialArgs = NULL, currentArgs = NULL) {
  if (is.null(specialArgs)) {
    stop("Please provide special arguments.")
  }
  
  necessaryCols <- list()
  currentCols <- names(currentArgs)
  
  if ('normalizedNodeCoverage' %in% names(specialArgs)) {
    if (!'run' %in% currentCols) {
      necessaryCols$run = NA
    }
    necessaryCols$node_coverage = NA
  }
  
  if ('biosampleTissue' %in% names(specialArgs)) {
    if (!'bio_sample' %in% currentCols) {
      necessaryCols$bio_sample = NA
    }
  }
  
  if ('STAT' %in% names(specialArgs)) {
    if (!'run' %in% currentCols) {
      necessaryCols$run = NA
    }
  }
  
  if ('librarySource' %in% names(specialArgs)) {
    if(!'run' %in% currentCols) {
      necessaryCols$run = NA
    }
  }
  
  return(necessaryCols)
}

retrieveParamInfo <- function(paramType = NA) {
  # Validate paramType
  validParamTypes <- c("categoricalPlotting", "numericalPlotting", 
                       "querying", "specialQuerying", "inequalityQuerying")
  if (!is.na(paramType) && !(paramType %in% validParamTypes)) {
    stop("Invalid paramType provided.")
  }
  
  # Return the list of parameters based on paramType
  switch(paramType,
         categoricalPlotting = getCategoricalPlottingParams(),
         numericalPlotting = getNumericalPlottingParams(),
         querying = getQueryingParams(),
         specialQuerying = getSpecialQueryingParams(),
         inequalityQuerying = getInequalityQueryingParams(),
         stop("A valid paramType must be specified."))
}

#' @title getAllDataProcessed
#'
#' @description Retrieves all data for the given query parameters processed by
#' Serratus (i.e. not just virus positive runs). 
#' 
#' @param conn A database connection object.
#' @param librarySource A character string specifying the library source taxa.
#' @param runID A vector of run IDs.
#' @param biosampleID A vector of biosample IDs.
#' @param bioprojectID A vector of bioproject IDs.
#' 
#' @return A tibble containing the processed data.
#' 
#' @examples
#' conn <- connectToDatabase()
#' query <- list(conn = conn, librarySource = "Nematoda")
#' data <- getAllDataProcessed(query)
getAllDataProcessed <- function(conn = NULL, librarySource = NA, runID = NA, biosampleID = NA, bioprojectID = NA) {
  if (is.null(conn) || !inherits(conn, "DBIConnection")) {
    stop("Invalid database connection.")
  }
  
  if (is.na(librarySource) && is.na(runID) && is.na(biosampleID) && is.na(bioprojectID)) {
    stop("Please provide at least one of librarySource, runID, biosampleID, or bioprojectID.")
  }
  
  librarySource <- if (!is.na(librarySource)) handleTaxSearch(librarySource) else NA
  queryFilters <- list(run = runID, bio_sample = biosampleID, bio_project = bioprojectID, scientific_name = librarySource)
  query <- constructQuery(filters = queryFilters, SQLtable = "srarun")
  
  tryCatch({
    data <- executeQuery(conn, query)
    data <- tibble(data)
    colNamesR <- SQLnames2Rnames(colnames(data))
    data <- setNames(data, colNamesR)
  }, error = function(e) {
    stop("Failed to execute query: ", e$message)
  })
  
  return(data)
}

#' @title getRunDiff 
#' 
#' @description Returns the rows in allData that are not present in virome.
#' Intended for use with the getAllDataProcessed function.
#' 
#' @param allData A data frame containing all data processed by Serratus for
#' a given query. See getAllDataProcessed. MUST contain a 'runID' column.
#' 
#' @param virome A data frame containing the virome data. MUST contain a 
#' 'runID' column.
#' 
#' @return A data frame containing the rows in allData that are not present in
#' virome.
#' 
#' @examples
#' conn <- connectToDatabase()
#' allData <- getAllDataProcessed(conn, runID = NA, librarySource = "Nematoda")
#' virome <- getViromeData(conn, runID = NA, librarySource = "Nematoda")
#' getRunDiff(allData, virome)
getRunDiff <- function(allData = NULL, virome = NULL) {
  if (is.null(allData) || is.null(virome) || !is.data.frame(allData) || !is.data.frame(virome)) {
    stop("allData and virome must be provided as data frames.")
  }
  
  requiredColumn <- "runID"
  if (!(requiredColumn %in% colnames(allData) && requiredColumn %in% colnames(virome))) {
    stop(sprintf("Both data frames must contain a '%s' column.", requiredColumn))
  }
  
  absentRuns <- setdiff(unique(allData[[requiredColumn]]), unique(virome[[requiredColumn]]))
  return(allData[allData[[requiredColumn]] %in% absentRuns, , drop = FALSE])
}

# [END]