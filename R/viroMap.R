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
viroMap <- function(virome = NULL, minCov = 2) {
  nodeCoverageMax <- virome %>%
    dplyr::group_by(sotu, run) %>%
    dplyr::filter(node_coverage > minCov) %>%
    dplyr::summarize(max_coverage = log(max(node_coverage, na.rm = TRUE)), .groups = "drop")

  # Get viral families
  families <- virome %>%
    dplyr::select(sotu, tax_family, node_coverage) %>%
    dplyr::filter(node_coverage > minCov) %>%
    dplyr::select(-node_coverage) %>%
    dplyr::distinct()

  # Get source species
  sourceSpecies <- virome %>%
    dplyr::select(run, scientific_name, node_coverage) %>%
    dplyr::filter(node_coverage > minCov) %>%
    dplyr::select(-node_coverage) %>%
    dplyr::distinct()

  # get row order
  rowOrder <- families %>%
    dplyr::arrange(desc(tax_family))
  rowOrder[is.na(rowOrder)] <- "Unknown"

  # get col order - drop species name
  sourceSpecies$scientific_name <- gsub(" .*",
                                        "", sourceSpecies$scientific_name)
  colOrder <- sourceSpecies %>%
    dplyr::arrange(desc(scientific_name))

  nodeCoverage <- nodeCoverageMax %>%
    pivot_wider(names_from = run, values_from = max_coverage) %>%
    replace(is.na(.), 0)

  # Sort
  nodeCoverage <- nodeCoverage %>%
    slice(match(rowOrder$sotu, nodeCoverage$sotu))
  nodeCoverage <- nodeCoverage[colOrder$run]

  colors <- colorRamp2(c(0, 10), c("white", "red"))
  plotMat <- as.matrix(nodeCoverage, byrow = TRUE)
  plotMat <- apply(plotMat, 2, as.numeric)

  rowAnnotation <- rowAnnotation(Virus_family=rowOrder$tax_family)
  colAnnotation <- HeatmapAnnotation(Source_species=colOrder$scientific_name)

  map <- Heatmap(plotMat,
                 name = "log(node_coverage)",
                 cluster_rows=TRUE, cluster_columns=TRUE,
                 show_row_names=FALSE, show_column_names=FALSE,
                 rect_gp = gpar(col = "white", lwd = 2),
                 width = unit(20, "cm"), height = unit(20, "cm"),
                 col = colors,
                 bottom_annotation=colAnnotation,
                 left_annotation=rowAnnotation)

  return(map)

}
