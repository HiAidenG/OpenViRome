#' @title viroMap
#' @description Plot a heatmap for the virome. Each column is a sample
#' and each row is an sOTU. The intensity of color for each cell is proportional
#' to the log of the node_coverage value for that sOTU in that sample. Samples
#' are grouped by source species. sOTUs are grouped by family.
#' @param virome A tbl_df containing the virome data.
#' @param minCov The minimum coverage value to include in the heatmap.
#' @importFrom dplyr group_by summarize filter distinct select
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly plot_ly layout
#' @export
#'
viroMap <- function(virome = NULL, minCov = 0) {

  virome <- virome[[1]]

  # Filter by minCov
  virome <- virome %>%
    dplyr::filter(node_coverage_norm > minCov)

  nodeCoverageMean <- virome %>%
    dplyr::group_by(sotu, run) %>%
    dplyr::summarize(mean_coverage = (max(node_coverage_norm, na.rm = TRUE)),
                     .groups = "drop")

  # Becomes the row names
  families <- virome %>%
    dplyr::select(sotu, tax_phylum) %>%
    dplyr::distinct()

  # Becomes the column names
  sourceSpecies <- virome %>%
    dplyr::select(run, scientific_name) %>%
    dplyr::distinct()

  # get row order
  rowOrder <- families %>%
    dplyr::arrange(desc(tax_phylum))
  rowOrder[is.na(rowOrder)] <- "Unknown"

  # get col order - drop species name
  sourceSpecies$scientific_name <- gsub(" .*",
                                        "", sourceSpecies$scientific_name)
  colOrder <- sourceSpecies %>%
    dplyr::arrange(desc(scientific_name))

  nodeCoverage <- nodeCoverageMean %>%
    pivot_wider(names_from = run, values_from = mean_coverage) %>%
    replace(is.na(.), 0)

  # Sort
  nodeCoverage <- nodeCoverage %>%
    slice(match(rowOrder$sotu, nodeCoverage$sotu))
  nodeCoverage <- nodeCoverage[colOrder$run]

  nodeCoverage <- log(nodeCoverage + 1)

  colors <- colorRamp2(c(0.0, max(nodeCoverage)), c("white", "#FF4040"))
  plotMat <- as.matrix(nodeCoverage, byrow = TRUE)
  plotMat <- apply(plotMat, 2, as.numeric)


  # Make row and column colors
  uniqueRows <- unique(rowOrder$tax_phylum)
  rowColors <- RColorBrewer::brewer.pal(length(uniqueRows), "Set3")
  names(rowColors) <- uniqueRows

  uniqueCols <- unique(colOrder$scientific_name)
  if (length(uniqueCols) < 3) {
    uniqueCols <- c(uniqueCols, uniqueCols, uniqueCols)
  }
  colColors <- RColorBrewer::brewer.pal(length(uniqueCols), "Dark2")
  names(colColors) <- uniqueCols


  rowAnnotation <- rowAnnotation(tax_phylum=rowOrder$tax_phylum,
                                    col = list(tax_phylum = rowColors))
  colAnnotation <- HeatmapAnnotation(Source_species=colOrder$scientific_name,
                                     col = list(Source_species = colColors))


  map <- Heatmap(plotMat,
                 name = "log(normalized node coverage)",
                 cluster_rows=TRUE, cluster_columns=TRUE,
                 show_row_names=FALSE, show_column_names=FALSE,
                 rect_gp = gpar(col = "white", lwd = 2),
                 width = unit(20, "cm"), height = unit(20, "cm"),
                 col = colors,
                 bottom_annotation=colAnnotation,
                 left_annotation=rowAnnotation
                 )
  return(map)

}
