# StatsModule.R

# my first idea for a statistical test I can do on this type of data is to
# compute the jaccard index for any two arbitrary 'viromes.' It would work as
# follows:
# - the user inputs two data frames of virome observations (which can be retrieved
# via the retrieveOTU function)
# - The function simply gathers all unique sOTUs observed in those viromes and
# considers the 'normalizedNodeCoverage' column as the observations for that sOTU,
# and then computes the jaccard index.

# Prerequisites: two viromes with >= 1 sotu, both containing the normalizedNodeCoverage column
getJaccard <- function(virome1 = NULL, virome2 = NULL) {
  # Check if the viromes are valid
  if (is.null(virome1) || is.null(virome2)) {
    stop("Please provide two viromes.")
  }

  # Check if the viromes are data frames and have at least one row
  if (!is.data.frame(virome1) || !is.data.frame(virome2)) {
    stop("Both viromes must be data frames.")
  }
  if (nrow(virome1) == 0 || nrow(virome2) == 0) {
    stop("Please provide viromes with at least one sOTU.")
  }
  if (!"normalizedNodeCoverage" %in% colnames(virome1) || !"normalizedNodeCoverage" %in% colnames(virome2)) {
    stop("Please provide viromes with the 'normalizedNodeCoverage' column.")
  }

  intersection <- getViromeIntersection(virome1 = virome1, virome2 = virome2, mode = 'sum')
  union <- getViromeUnion(virome1 = virome1, virome2 = virome2, mode = 'sum')
  jaccardIndex <- intersection / union

  return(jaccardIndex)
}


getViromeIntersection <- function(virome1 = NULL, virome2 = NULL, mode = 'sum') {
    if (mode != 'sum' && mode != 'sotu') {
        stop("Please provide a valid mode: 'sum' or 'sotu'.")
    }

    virome1 <- collapseSotuObservations(virome1)
    virome2 <- collapseSotuObservations(virome2)

    sotus <- unique(c(virome1$sOTU, virome2$sOTU))

    # Consider the 'normalizedNodeCoverage' column as the observations for that sOTU
    virome1Obs <- virome1$normalizedNodeCoverage[match(sotus, virome1$sOTU)]
    virome2Obs <- virome2$normalizedNodeCoverage[match(sotus, virome2$sOTU)]

    # Replace NAs with zeros
    virome1Obs[is.na(virome1Obs)] <- 0
    virome2Obs[is.na(virome2Obs)] <- 0

    if (mode == 'sum') {
        intersection <- sum(pmin(virome1Obs, virome2Obs))
    } else if (mode == 'sotu') {
        # Create a data frame of all overlapping sOTUs, their normalized node coverage in virome1 and virome2
        intersection <- tibble(
            sotu = sotus,
            virome1Coverage = virome1Obs,
            virome2Coverage = virome2Obs,
        )
        # Keep only rows with non-zero values in both coverages to represent overlapping sOTUs
        intersection <- intersection[intersection$virome1Coverage > 0 & intersection$virome2Coverage > 0, ]
        # Make a normalizedNodeCoverage column that's the sum of the two coverages
        intersection$normalizedNodeCoverage <- intersection$virome1Coverage + intersection$virome2Coverage
    }

    return(intersection)
}



getViromeUnion <- function(virome1 = NULL, virome2=NULL, mode = 'sum') {
    if (mode != 'sum' && mode != 'sotu') {
        stop("Please provide a valid mode: 'sum' or 'sotu'.")
    }
    # Get unique sotus. If there are multiple observations of the same sotu, take the sum
    virome1 <- collapseSotuObservations(virome1)
    virome2 <- collapseSotuObservations(virome2)

    sotus <- unique(c(virome1$sOTU, virome2$sOTU))

    # Consider the 'normalizedNodeCoverage' column as the observations for that sOTU
    virome1Obs <- virome1$normalizedNodeCoverage[match(sotus, virome1$sOTU)]
    virome2Obs <- virome2$normalizedNodeCoverage[match(sotus, virome2$sOTU)]

    # Replace NAs with zeros
    virome1Obs[is.na(virome1Obs)] <- 0
    virome2Obs[is.na(virome2Obs)] <- 0

    if (mode == 'sum') {
        union <- sum(pmax(virome1Obs, virome2Obs))
    } else if (mode == 'sotu') {
        # return a data frame of all sOTUs and their normalizedNodeCoverage
        union <- tibble(sotu = sotus, normalizedNodeCoverage = pmax(virome1Obs, virome2Obs))
    }

    return(union)
}

# Next test (?)

# a virome network is a graph of viromes where the nodes are viromes and the edges
# are the jaccard index between the viromes. This function will take a list of viromes
# and a similarity cutoff and return a dataframe of the edges of the graph.
# viromeList is a named list of tibbles
makeViromeNetwork <- function(viromeList = NULL, similarityCutoff = NULL) {
  if (is.null(viromeList)) {
    stop("Please provide a list of viromes.")
  }

  if (is.null(similarityCutoff)) {
    stop("Please provide a similarity cutoff.")
  }
  if (similarityCutoff < 0 || similarityCutoff > 1) {
    stop("Please provide a similarity cutoff between 0 and 1.")
  }

  # Initialize an empty dataframe to store the network edges
  networkEdges <- tibble(virome1 = character(), virome2 = character(), jaccardIndex = numeric())

  # Compute the jaccard index for all pairs of viromes
  for (i in 1:(length(viromeList) - 1)) {
    for (j in (i + 1):length(viromeList)) {
      virome1 <- viromeList[[i]]
      virome2 <- viromeList[[j]]
      jaccardIndex <- getJaccard(virome1 = virome1, virome2 = virome2)

      if (jaccardIndex >= similarityCutoff) {
        networkEdges <- rbind(networkEdges, tibble(virome1 = names(viromeList)[i], virome2 = names(viromeList)[j], jaccardIndex = jaccardIndex))
      }
    }
  }

  # Filter the edges based on the similarity cutoff
  networkEdges <- networkEdges[networkEdges$jaccardIndex >= similarityCutoff, ]

  return(networkEdges)
}

collapseSotuObservations <- function(virome = NULL) {
  virome <- virome %>% group_by(sOTU) %>% summarise(normalizedNodeCoverage = sum(normalizedNodeCoverage))

  return(virome)
}

addZScoreColumn <- function(data) {
  # Calculate the mean and standard deviation of normalizedNodeCoverage
  mean_coverage <- mean(data$normalizedNodeCoverage)
  sd_coverage <- sd(data$normalizedNodeCoverage)

  # Calculate the z-score for each sotu
  z_scores <- (data$normalizedNodeCoverage - mean_coverage) / sd_coverage

  # Add the z-scores as a new column to the data
  data$zScore <- z_scores

  return(data)
}

getViromeSummary <- function(virome = NULL, serratusProcessed = NULL) {
  # Get the number of unique sOTUs in the virome
  sOTUs <- virome %>% distinct(sOTU) %>% nrow()
  
  # Get the total number of runs processed 
  runsProcessed <- serratusProcessed %>% distinct(runID) %>% nrow()
  
  # Get the number of virus positive runs
  virusPositiveRuns <- virome %>% distinct(runID) %>% nrow()
  
  # Mean normalized node coverage
  meanCoverage <- mean(virome$normalizedNodeCoverage)
  
  # Max normalized node coverage
  maxCoverage <- max(virome$normalizedNodeCoverage)
  
  # Number of distinct tissue types
  tissueTypes <- virome %>% distinct(biosampleTissue) %>% nrow()
  
  # Avg genBank identity
  avgGBpID <- mean(virome$GBpID)
  
  # Number of distinct taxonomic sources
  sourceSpecies <- virome %>% distinct(librarySource) %>% nrow()
  
  bioprojects <- virome %>% distinct(bioprojectID) %>% nrow()
  
  # Rounding numeric values to the nearest hundredth
  meanCoverage <- round(meanCoverage, 2)
  maxCoverage <- round(maxCoverage, 2)
  avgGBpID <- round(avgGBpID, 2)
  
  return(c(sOTUs = sOTUs, runsProcessed = runsProcessed, 
           virusPositiveRuns = virusPositiveRuns, 
           meanCoverage = meanCoverage, 
           maxCoverage = maxCoverage, 
           tissueTypes = tissueTypes, 
           avgGBpID = avgGBpID,
           bioprojects = bioprojects,
           sourceSpecies = sourceSpecies))
}



#
