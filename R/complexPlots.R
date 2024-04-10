# complexPlots.R
# This module provides several functions for creating 'complex' i.e. non-standard plots for virome data.
# Many of these have stringent preconditions on the data that must be met before they can be plotted.

# Plot 1: Barplot of the number of virus positive SRA runs over all runs processed per library source
# Prereq: virome object must contain the columns 'runID' and 'librarySource'
# Takes in a virome data frame and a an allData data frame (see getAllDataProcessed())
# Plots a stacked barchart with the proportion of virus positive runs per library source
# Uses the getRunDiff() function to calculate the number of virus positive runs per library source
plotVirusPositive <- function(viromeData = NULL, processedRunData = NULL, title = 'Virus positive runs by library Source') {
	if (is.null(viromeData)) {
		stop("No virome data provided.")
	}

	if (is.null(processedRunData)) {
		stop("No processed runs data provided.")
	}

	if (!"runID" %in% colnames(viromeData)) {
		stop("The virome object must contain a 'runID' column.")
	}

	if (!"librarySource" %in% colnames(viromeData)) {
		stop("The virome object must contain a 'librarySource' column.")
	}

	# Create a data frame with columns for runID, librarySource, and virusPositive status
	run_data <- tibble::tibble(
		runID = processedRunData$runID,
		librarySource = processedRunData$librarySource,
		virusPositive = ifelse(processedRunData$runID %in% viromeData$runID, "Virus Positive", "Virus Negative")
	)

	# Create the ggplot object for a horizontal bar chart
	ggplot_object <- ggplot2::ggplot(run_data, aes(x = librarySource, fill = virusPositive)) +
		geom_bar(stat = "count") +
		scale_fill_manual(values = c("Virus Positive" = "#CC79A7", "Virus Negative" = "#999999"), guide = FALSE) +
		labs(title = title, x = "Library Source", y = "Number of Runs") +
		coord_flip() +
		theme_bw() +
		theme_cowplot(12)

	# Convert the ggplot object to a Plotly object
	plotly_object <- ggplotly(ggplot_object)
	
	return(plotly_object)
}


# Plot 2: sankey diagram for flow of GB acc to sotu to run to bio
# prereq: virome must contain the columns: 'runID', 'sOTU', 'GBAcc,' 'bioprojectID'
# Takes in a virome data frame with above columns, plots a sankey diagram
plotSankey <- function(virome, linkColors = NULL) {
  if (is.null(virome)) {
    stop("No virome data provided.")
  }
  
  preReqCols <- c("sOTU", "GBAcc", "bioprojectID", 'librarySource')
  if (!all(preReqCols %in% colnames(virome))) {
    stop("The virome object must contain columns: 'sOTU', 'GBAcc', 'bioprojectID', and 'librarySource'")
  }
  
  virome <- virome %>%
    dplyr::select(sOTU, bioprojectID, GBAcc, librarySource)
  
  # Assign colors to bioprojectID, use provided linkColors if available
  if (is.null(linkColors)) {
    bioColors <- setNames(grDevices::rainbow(length(unique(virome$bioprojectID))), unique(virome$bioprojectID))
  } else {
    if (length(linkColors) != length(unique(virome$bioprojectID))) {
      stop("The length of linkColors should match the number of unique bioprojectIDs.")
    }
    bioColors <- setNames(linkColors, unique(virome$bioprojectID))
  }
  
  virome$color <- bioColors[virome$bioprojectID]
  
  # Create the links for the Sankey diagram
  links <- virome %>% 
    dplyr::group_by(GBAcc, sOTU, bioprojectID) %>%
    dplyr::summarise(value = n(), color = first(color)) %>%
    dplyr::rename(source = GBAcc, target = sOTU) %>%
    dplyr::bind_rows(virome %>%
                       dplyr::group_by(sOTU, bioprojectID) %>%
                       dplyr::summarise(value = n(), color = first(color)) %>%
                       dplyr::rename(source = sOTU, target = bioprojectID)) %>%
    dplyr::bind_rows(virome %>%
                       dplyr::group_by(bioprojectID, librarySource) %>%
                       dplyr::summarise(value = n(), color = first(color)) %>%
                       dplyr::rename(source = bioprojectID, target = librarySource))
  
  # Create a list of unique nodes
  nodes <- data.frame(name = unique(c(links$source, links$target)))
  
  # Convert node names to indices
  links$source <- match(links$source, nodes$name) - 1
  links$target <- match(links$target, nodes$name) - 1
  
  p <- plot_ly(
    type = 'sankey',
    domain = list(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    node = list(
      pad = 10,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      ),
      label = nodes$name,
      color = rep("grey", length(nodes$name))  # Node Coloring
    ),
    link = list(
      source = links$source,
      target = links$target,
      value = links$value,
      color = links$color  # Link Coloring by bioprojectID or user-defined colors
    )
  )
  
  p  # return the plotly object
}

# Plot 3: Virome network plot
# Break down the virome into individual bioprojects and their associated sOTUs.
# then, create a graph network where each node is a bioproject.
# edges are drawn between bioprojects based on a jacccard similarity of sOTUs (makeViromeNetwork)
# Prereq: virome must contain the columns: 'sOTU', 'bioprojectID' and ideally consist of
# multiple bioprojects
#' @title plotViromeNetwork
#' 
#' @description
#' 
#' @param virome A data frame containing columns 'sOTU,' 'bioprojectID,' and 
#' 'normalizedNodeCoverage'
#' 
#' @param threshold A numeric value between 0 and 1 that determines the threshold for
#' the Jaccard similarity of sOTUs between bioprojects. Default is 0.0.
#' 
#' @param splitVariable A character string that determines how the virome data should be
#' split. Default is 'bioprojectID.'
#' 
#' @return A networkD3 plot
#' 
#' @examples
#' virome <- retrieveOTUs(bioprojectID = NA, librarySource = "Tylenchoidea", 
#' normalizedNodeCoverage = NA)
#' plotViromeNetwork(virome = virome, threshold = 0.01, splitVariable = "bioprojectID")
#' 
#' @importFrom networkD3 simpleNetwork
plotViromeNetwork <- function(virome = NULL, threshold = 0.01, splitVariable = "bioprojectID") {
  if (is.null(virome)) {
    stop("No virome data provided.")
  }
  
  allowedSplits <- c("bioprojectID", "librarySource")
  
  preReqCols <- c("sOTU", splitVariable)
  if (!all(preReqCols %in% colnames(virome))) {
    stop("The virome object must contain columns: 'sOTU' and '", splitVariable, "'.")
  }
  
  # Split the virome into a list of bioprojects and their associated sOTUs
  viromeList <- split(virome, virome[[splitVariable]])
  
  # Make a call to makeViromeNetwork
  network <- makeViromeNetwork(viromeList, threshold)
  
  # Create the network plot
  networkPlot <- networkD3::simpleNetwork(network)
  
  return(networkPlot)
}

validateSet <- function(set) {
  if (!is.character(set) || nchar(set) == 0) {
    stop("Each set must be a non-empty string.")
  }
}

filterViromeData <- function(virome, set) {
  if (!is.data.frame(virome)) {
    stop("virome must be a data frame.")
  }
  if (!any(c("bioprojectID", "librarySource") %in% colnames(virome))) {
    stop("The virome object must contain 'bioprojectID' or 'librarySource' columns.")
  }
  filtered <- virome[virome$bioprojectID == set | virome$librarySource == set, ]
  if (nrow(filtered) == 0) {
    stop("No matching entries found in virome for the provided set.")
  }
  filtered
}

plotsOTUVenn <- function(virome = NULL, set1 = NULL, set2 = NULL, set3 = NA, set4 = NA) {
  if (is.null(virome)) {
    stop("No virome data provided.")
  }
  if (is.null(set1) || is.null(set2)) {
    stop("At least two sets must be provided.")
  }
  
  # Prepare the list of sets
  sets <- c(set1, set2, set3, set4)
  setNames <- c(set1, set2, set3, set4)
  validSets <- sets[!is.na(sets)]
  setNames <- setNames[!is.na(setNames)]
  
  lapply(validSets, validateSet)
  
  # Filter and process the virome data
  sOTUs <- setNames(validSets, setNames) # Assign names to the sets
  filteredSets <- lapply(sOTUs, function(set) filterViromeData(virome, set))
  
  # Prepare the sOTUs for the Venn diagram
  sOTUsData <- lapply(filteredSets, function(set) {
    if (!is.null(set) && "sOTU" %in% colnames(set)) {
      as.integer(sub("^u", "", set$sOTU))
    }
  })
  
  if (length(sOTUsData) < 2 || length(sOTUsData) > 4) {
    stop("The number of sets for the Venn diagram must be between 2 and 4.")
  }
  
  # Create the Venn diagram with proper set names
  ggvenn::ggvenn(sOTUsData, names(sOTUsData), text_size = 4, set_name_size =4)
}

# Plot 4: PCA of virome diversity
# Step 1: Split the virome on either a) bioprojectID or b) librarySource
# Step 2: Work out sOTU diversity for each bioproject/librarySource. i.e. for
# every bioproject/librarySource, make a dataframe with sOTUs as rows and
# bioproject/librarySource as columns. The values in the dataframe are the
# normalizedNodeCoverage of the sOTUs. This will be a sparse matrix.
# Step 4: PERFORM PCA on the sparse matrix
# Step 5: Plot the PCA
plotViromePCA <- function(virome = NULL, mode = 'bioprojectID', sOTU = FALSE, title = '') {
  if (is.null(virome) || !is.data.frame(virome)) {
    stop("virome must be a non-empty data frame.")
  }
  
  if (!(mode %in% c('bioprojectID', 'librarySource'))) {
    stop("mode must be either 'bioprojectID' or 'librarySource'.")
  }
  
  if (!mode %in% names(virome)) {
    stop(paste("The column", mode, "does not exist in the virome data frame."))
  }
  
  # Aggregate normalizedNodeCoverage by sOTU and mode
  aggregated_data <- virome %>%
    group_by(across(all_of(mode)), sOTU) %>%
    summarise(normalizedNodeCoverage = sum(normalizedNodeCoverage), .groups = 'drop')
  
  # Convert to wide format
  diversity_data <- aggregated_data %>%
    pivot_wider(names_from = all_of(mode), values_from = normalizedNodeCoverage, values_fill = list(normalizedNodeCoverage = 0))

  if (!sOTU) {
    diversity_data <- t(diversity_data)
    diversity_data <- diversity_data[-1,]
    annotation <- rownames(diversity_data)
  }
  else {
    annotation <- diversity_data$sOTU
    diversity_data <- diversity_data[-1]
  }

  # Convert to numeric
  diversity_data <- as.data.frame(apply(diversity_data, 2, as.numeric))
  # Perform PCA
  pca_results <- prcomp(as.matrix(diversity_data), center = TRUE, scale. = TRUE)
  
  # Prepare the data for plotting
  pca_data <- as.data.frame(pca_results$x)
  pca_data$annotation <- annotation

  # Create the interactive plot with proper identifiers
  pca_plot <- plot_ly(pca_data, x = ~PC1, y = ~PC2, text = ~annotation, hoverinfo = 'text+x+y', type = 'scatter', mode = 'markers') %>%
    layout(title = title,
           xaxis = list(title = "Principal Component 1"),
           yaxis = list(title = "Principal Component 2"),
           hovermode = 'closest')
  
  return(pca_plot)
}
