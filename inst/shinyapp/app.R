# app.R

library(shiny)
library(plotly)
library(tidyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(bs4Dash)
require(circlize)
require(grDevices)
require(fresh)
require(cowplot)
library(tidyr)

viromeTheme <- fresh::create_theme(
  bs4dash_status(
    primary = "#6495ED",
    secondary = "#818589",
    success = "#00CD66",
  )
)

ui <- bs4DashPage(
  freshTheme = viromeTheme,
  title = "openviRome",
  header = bs4Dash::dashboardHeader(
    title = "openviRome v0.02",
    titleWidth = 250,
    rightUi = uiOutput("downloadButtonUI")  # Placeholder for download button
  ),
  sidebar = bs4DashSidebar(
    skin = "dark",
    sidebarMenu(
      id = "tabs",
      menuItem(
        text = "Home",
        tabName = "homeTab",
        icon = icon("home")
      ),
      menuItem(
        text = "Overview",
        tabName = "ViromeOverview",
        icon = icon("chart-pie")
      ),
      menuItem(
        text = "Complex Plots",
        tabName = "ComplexPlots",
        icon = icon("chart-bar")
      ),
      # menuItem(
      #   text = "Source-sOTU Heatmap",
      #   tabName = "sOTUHeatmap",
      #   icon = icon("vial-virus")
      # ),
      menuItem(
        text = "Virome Network",
        tabName = "vNetwork",
        icon = icon("viruses")
      )
    )
  ),
  body = bs4DashBody(
    bs4Dash::bs4TabItems(
      bs4Dash::bs4TabItem(
        tabName = "homeTab",
        fluidRow(
          column(
            width = 6, offset = 0,
            h2("openviRome"),
            p("Explore Serratus data through a viromics lens! openviRome aims to provide a user-friendly interface for virus discovery by democratizing public sequencing data.")
          )
        ), # Row 1 close
        fluidRow(
          bs4Dash::box(title = "Query Serratus",
                       p("Use the form below to construct your query to Serratus. For more information about Serratus queries, see the 'help' tab."),
                       closable = FALSE,
                       collapsible = FALSE,
                       width = 6,
                       solidHeader = TRUE,
                       status = "primary",
                       # Input fields
                       p("Find OTUs appearing in..."),
                       textInput("SRARuns", label = NULL,
                                 placeholder = "SRA Accessions (comma separated)",
                                 width = '100%'),
                       textInput("BioSamples", label = NULL,
                                 placeholder = "BioSample Accessions (comma separated)",
                                 width = '100%'),
                       textInput("BioProjects", label = NULL,
                                 placeholder = "BioProject Accessions (comma separated)",
                                 width = '100%'),
                       p("With the following metadata annotations..."),
                       textInput("sourceOrganism", label = NULL,
                                 placeholder = "Source Organism (taxonomic name)",
                                 width = '50%'),
                       textInput("tissueType", label = NULL,
                                 placeholder = "Tissue Type",
                                 width = '50%'),
                       p("Then filter for observations with..."),
                       textInput("minCoverage", label = NULL,
                                 placeholder = "Minimum Coverage",
                                 width = '50%'),
                       # textInput("STATTaxa", label = NULL,
                       #           placeholder = "STAT Taxa (comma separated)",
                       #           width = '50%'),
                       # textInput("minKmerPerc", label = NULL,
                       #           placeholder = "min kmer perc",
                       #           width = '50%'),
                       p("And only show sOTUs with the following predicted taxonomy..."),
                       textInput("Phylum", label = NULL,
                                 placeholder = "Phylum",
                                 width = '50%'),
                       textInput("Order", label = NULL,
                                 placeholder = "Order",
                                 width = '50%'),
                       textInput("Family", label = NULL,
                                 placeholder = "Family",
                                 width = '50%'),
                        textInput("Genus", label = NULL,
                                 placeholder = "Genus",
                                 width = '50%'),
                        textInput("Species", label = NULL,
                                  placeholder = "Species",
                                  width = '50%'),
                       actionButton("submitQuery", "Submit", style = "width: 100%;")
          ),
          bs4Dash::box(title = "Upload Virome",
                       p("Upload a formatted virome file in CSV format."),
                       p("See the 'help' tab for more information."),
                       p("Example .csv provided at openviRome/inst/extdata/Meloidogyne_virome.csv"),
                       closable = FALSE,
                       collapsible = FALSE,
                       width = 6,
                       solidHeader = FALSE,
                       status = "primary",
                       fileInput("uploadFile", label = NULL,
                                 placeholder = "Upload Virome",
                                 width = '50%'),
                       actionButton("submitFile", "Submit", style = "width: 50%;")
          )
        ), # Row 2 close
        fluidRow(
          column(width = 12,
                 # Display loading animation
                 uiOutput("statsCard")
          )
        ), # Row 3 close
        fluidRow(
          bs4Dash::bs4Card(
            title = "Help",
            width = 12,
            tabPanel(
            p("This tool is currently still in beta. Please report any bugs or issues to the openviRome GitHub page."),
            h4("How to use openviRome"),
            p("openviRome is designed to be able to provide a 'virome-style' analysis for any arbitrary set of SRA data. As of v0.02, openviRome can be used in three ways:"),
            p("1. Query Serratus by taxon. This will construct a virome from all SRA data processed by Serratus that are annotated as belonging to the taxon specified. This approach has limitations, notably in that SRA metadata annotations are not particularly reliable. Please also note that this relies on the NCBI taxonomy, so taxonomic names must be NCBI-recognized. If you are unsure, check the ", a("NCBI Taxonomy Browser", href = "https://www.ncbi.nlm.nih.gov/taxonomy"), "."),
            p("2. Query Serratus by SRA accession. The virome will consist of all SRA accessions provided in your comma-separated list if and only if they have been processed by Serratus. At the time of writing (December 8th, 2023) Serratus v2 has processed X runs (xx%)."),
            p("3. Upload a virome file. This file should be a .csv with the following columns:"),
            p("    - run: SRA accession"),
            p("    - scientific_name: source organism taxonomic name"),
            p("    - bio_sample: SRA biosample accession"),
            p("    - sotu: source-sOTU ID"),
            p("    - gb_acc: genbank accession of closest alignemnt."),
            p("    - gb_pid: percent identity of closest alignment."),
            p("    - tax_species: species name of closest alignment."),
            p("    - tax_phylum: phylum name of closest alignment."),
            p("    - node_coverage_norm: read coverage for sOTU normalized by library size."),
            p("For more information about Serratus, see the ", a("wiki", href = "https://github.com/ababaian/serratus/wiki"), ".")
            )
          )
        ) # Row 4 close
      ), # Close homeTab
      bs4TabItem(
        tabName = "ViromeOverview",
        h2("Virome Overview"),
        fluidRow(
          column(
            width = 12,
            uiOutput("dynamicSlider")
          )
        ),
        fluidRow(
          column(
            width = 6,
            bs4Dash::box(
              title = "Histogram",
              status = "primary",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = NULL,
              plotlyOutput("HistPlot"),
              selectInput("histVariableSelect", 
                          "Choose a Variable", 
                          choices = retrieveParamInfo(paramType = "numericalPlotting"), 
                          selected = "normalizedNodeCoverage" 
              )
            )
          ),
          column(
            width = 6,
            bs4Dash::box(
              title = "sOTU Predicted Taxonomy",
              status = "primary",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = NULL,
              plotlyOutput("PieChartPlot"),
              selectInput("pieVariableSelect", 
                          "Choose a Variable", 
                          choices = c("phylum", "order", "family", "genus", "species"), 
                          selected = "phylum" 
              )
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            bs4Dash::box(
              title = "Violin Plot",
              status = "primary",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = NULL,
              plotlyOutput("ViolinPlot"),
              selectInput("violinXVariableSelect", 
                          "Choose a Variable for X-axis", 
                          choices = retrieveParamInfo(paramType = "numericalPlotting"), 
                          selected = "normalizedNodeCoverage" 
              ),
              selectInput("violinYVariableSelect", 
                          "Choose a Variable for Y-axis", 
                          choices = retrieveParamInfo(paramType = "categoricalPlotting"), 
                          selected = "phylum" 
              )
            )
          )
        ),
        # fluidRow(
        #   column(
        #     width = 12,
        #     bs4Dash::box(
        #       title = "Scatter Plot",
        #       status = "primary",
        #       solidHeader = FALSE,
        #       collapsible = FALSE,
        #       width = NULL,
        #       plotlyOutput("ScatterPlot"),
        #       selectInput("scatterXVariableSelect", 
        #                   "Choose a Variable for X-axis", 
        #                   choices = retrieveParamInfo(paramType = "numericalPlotting"), 
        #                   selected = "nodeCoverageNorm" 
        #       ),
        #       selectInput("scatterYVariableSelect", 
        #                   "Choose a Variable for Y-axis", 
        #                   choices = retrieveParamInfo(paramType = "numericalPlotting"), 
        #                   selected = "GBpID" 
        #       )
        #     )
        #   )
        # ),
        fluidRow(
          column(
            width = 12,
            bs4Dash::box(
              title = "Sankey Plot",
              status = "primary",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = NULL,
              plotlyOutput("sankeyPlot")
            )
          )
        )
      ), # Close ViromeOverview
      bs4TabItem(
        tabName = "ComplexPlots",
        h2("Additional Plots"),
        fluidRow(
          column(
            width = 12,
            bs4Dash::box(
              title = "Virus Positivity",
              plotlyOutput("plotVirusPositive"),
              collapsible = FALSE,
              width = NULL
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            bs4Dash::box(
              title = "PCA",
              plotlyOutput("plotPCA"),
              collapsible = FALSE,
              width = NULL,
              selectInput("pcaVariable",
                          "Choose a Variable",
                          choices = c("bioprojectID", "librarySource"),
                          selected = "bioprojectID"
              ),
              selectInput("pcaType",
                          "Choose whether to plot sOTUs (TRUE) or samples (FALSE)",
                          choices = c('TRUE', 'FALSE'),
                          selected = 'FALSE')
            
          )
        )
       )
      ), # Close complex plots
      bs4TabItem(
        tabName = "vNetwork",
        h2("sOTU Graph Network"),
        fluidRow(
          column(
            width = 12,
            bs4Dash::box(
              title = "sOTU Network",
              networkD3::simpleNetworkOutput("viromeNetwork"),
              collapsible = FALSE,
              width = NULL,
              selectInput("networkSplitVariable",
                          "Choose the nodes",
                          choices = c("bioprojectID", "librarySource"),
                          selected = "bioprojectID"),
              sliderInput("networkCutoff",
                          "Choose the cutoff for Jaccard similarity",
                          min = 0,
                          max = 1,
                          value = 0.1,
                          step = 0.01)
              )
            )
          )
        )
    )
  ), # Close body
  controlbar = NULL,
  footer = bs4DashFooter(),
  help = NULL,
  dark = NULL # Disable these to clean up the UI
)

server <- function(input, output, session) {

  viromeData <- reactiveVal(NULL)
  serratusProcessed <- reactiveVal(NULL)
  selectedPhylum <- reactiveVal(NULL)
  
  # ======================= Query button ============================
  observeEvent(input$submitQuery, {
    showModal(modalDialog(
      title = "Processing...",
      "Your data is being loaded. Depending on the size of your query, this may take a few minutes.",
      "You can click anywhere outside of this dialog to close it.",
      easyClose = TRUE,
      footer = NULL
    ))

    errorOccurred <- FALSE
    
    # Connect to the server
    conn <- connectToDatabase()
    
    # Check that the connection was successful
    if (is.null(conn)) {
      errorOccurred <- TRUE
      showModal(modalDialog(
        title = "Error",
        "An error occurred while connecting to the database. Please try again.",
        easyClose = TRUE,
        footer = NULL
      ))
    }

    # Check each of the input fields for data
    query <- list()
    if (input$SRARuns != "") {
      query$runID <- input$SRARuns
    }
    # Otherwise, make runID NA
    else {
      query$runID <- NA
    }

    if (input$BioSamples != "") {
      query$biosampleID <- input$BioSamples
    } else {
      query$biosampleID <- NA
    }

    if (input$BioProjects != "") {
      query$bioprojectID <- input$BioProjects
    } else {
      query$bioprojectID <- NA
    }

    if (input$sourceOrganism != "") {
      query$librarySource <- input$sourceOrganism
    } else {
      query$librarySource <- NA
    }

    if (input$tissueType != "") {
      query$biosampleTissue <- input$tissueType
    } else {
      query$biosampleTissue <- NA
    }

    if (input$minCoverage != "") {
      query$normalizedNodeCoverage <- as.numeric(input$minCoverage)
    } else {
      query$normalizedNodeCoverage <- 0
    }

    # TODO: Fix this! STAT searches are currently broken
    # if (input$STATTaxa != "") {
    #   query$STAT <- input$STATTaxa
    # } else {
    #   query$STAT <- NA
    # }

    # if (input$minKmerPerc != "") {
    #   query$kmerPerc <- as.numeric(input$minKmerPerc)
    # } else {
    #   query$kmerPerc <- 0
    # }

    if (input$Phylum != "") {
      query$phylum <- input$Phylum
    } else {
      query$phylum <- NA
    }

    if (input$Order != "") {
      query$order <- input$Order
    } else {
      query$order <- NA
    }

    if (input$Family != "") {
      query$family <- input$Family
    } else {
      query$family <- NA
    }

    if (input$Genus != "") {
      query$genus <- input$Genus
    } else {
      query$genus <- NA
    }

    if (input$Species != "") {
      query$species <- input$Species
    } else {
      query$species <- NA
    }
    
    query$GBpID <- NA
    query$GBAcc <- NA

    # Retrieve the virome data
    args <- c(conn = conn, query)
    virome <- do.call("retrieveOTUs", args)
    virome <- addZScoreColumn(virome)
    
    serratusProcessed <- getAllDataProcessed(conn = conn, librarySource = query$librarySource,
                                             runID = query$runID, biosampleID = query$biosampleID,
                                             bioprojectID = query$bioprojectID)
    
    if (!errorOccurred) {
      viromeData(virome)
      serratusProcessed(serratusProcessed)
      removeModal()
    }
  })

  # TODO: Update this. 
  # ======================= Upload button ============================
  observeEvent(input$submitFile, {

    req(input$uploadFile)
    virome <- read.csv(input$uploadFile$datapath)
    virome <- tibble::as_tibble(virome)

    if (viromeFormatCheck(virome = virome)) {
      showModal(modalDialog(
        title = "Processing...",
        "Retrieving run data. This should only take a moment.",
        easyClose = TRUE,
        footer = NULL
      ))

      removeModal()
    } else {
      showNotification("Error: Your file is not in the correct format.", type = "error")
    }

  })

  # ======================= Download button ============================
  output$downloadButtonUI <- renderUI({
    req(viromeData())
    downloadButton("downloadVirome", "Download Virome CSV")
  })

  # Define the download handler
  output$downloadVirome <- downloadHandler(
    content = function(file) {
      data <- viromeData()
      write.csv(data, file, row.names = FALSE)
    },
    filename = function() {
      paste("virome_", Sys.Date(), ".csv", sep = "")
    }
  )
  
  # ======================= Summary stats ============================
  output$statsCard <- renderUI({

    if (!is.null(viromeData()) && !is.null(serratusProcessed())) {
      stats <- getViromeSummary(virome = viromeData(), serratusProcessed = serratusProcessed())
      # Create a bs4Card to display the stats
      bs4Card(
        title = "Virome Summary Stats",
        status = "success",
        solidHeader = TRUE,
        collapsible = FALSE,
        width = NULL,
        # Display the stats in the card
        div(
          style = "display: flex; flex-wrap: wrap; text-align: center;",
          div(style = "flex: 33%; padding: 5px;", p("Number of Runs Processed: ", stats["runsProcessed"])),
          div(style = "flex: 33%; padding: 5px;", p("Number of Virus Positive Runs: ", stats["virusPositiveRuns"])),
          div(style = "flex: 33%; padding: 5px;", p("Distinct sOTUs: ", stats["sOTUs"])),
          div(style = "flex: 33%; padding: 5px;", p("Mean Read Coverage: ", stats["meanCoverage"])),
          div(style = "flex: 33%; padding: 5px;", p("Max Read Coverage: ", stats["maxCoverage"])),
          div(style = "flex: 33%; padding: 5px;", p("Distinct Tissue Types: ", stats["tissueTypes"])),
          div(style = "flex: 33%; padding: 5px;", p("Distinct Bioprojects: ", stats["bioprojects"])),
          div(style = "flex: 33%; padding: 5px;", p("Distinct Source Organisms: ", stats["sourceSpecies"])),
          div(style = "flex: 33%; padding: 5px;", p("Mean GenBank Identity: ", stats["avgGBpID"]))
        )
      )
    } else {
      #Display a greyed-out box with "No data loaded" message
      bs4Card(
        title = "Virome Summary Stats",
        status = "secondary",
        solidHeader = TRUE,
        collapsible = FALSE,
        width = NULL,
        # Display text in center of card
        div(
          style = "text-align: center;",
          p("No data loaded")
        )
      )
    }
  })
  
  output$dynamicSlider <- renderUI({
    req(viromeData())
    virome <- viromeData()
    
    sliderInput("zScoreSlider", 
                "Filter by sOTU coverage z-score:", 
                min = floor(min(virome$zScore, na.rm = TRUE)), 
                max = ceiling(max(virome$zScore, na.rm = TRUE)), 
                value = min(virome$zScore, na.rm = TRUE),
                step = 0.1)
  })
  
  # Create a reactive expression for the filtered data
  plotData <- reactive({
    req(viromeData())
    virome <- viromeData()
    minZScore <- input$zScoreSlider
    
    # Check if the slider is at its minimum value, indicating no filtering
    if (minZScore == min(virome$zScore, na.rm = TRUE)) {
      return(virome)  # Return the unfiltered data
    } else {
      # Otherwise, filter the data based on the selected minimum zScore
      return(virome[virome$zScore >= minZScore, ])
    }
  })

  # ======================= Pie Chart ============================
  output$PieChartPlot <- renderPlotly({
    req(plotData())  # Ensure data is available
    
    pieSelectedVariable <- input$pieVariableSelect  # Get the selected variable from the input
    req(pieSelectedVariable)  # Ensure a variable is selected
    
    filtered <- plotData()
    # Get unique sOTUs
    filtered <- filtered %>% distinct(sOTU, .keep_all = TRUE)
      
    
    plotPie(virome = filtered, variable = pieSelectedVariable, title = NA)  # Use the selected variable for plotting
  })
  
  # ======================= Histogram ============================
  output$HistPlot <- renderPlotly({
    req(plotData())  # Ensure data is available
    
    histSelectedVariable <- input$histVariableSelect  
    req(histSelectedVariable) 
    
    plotData <- plotData()
    filtered <- plotData %>% distinct(sOTU, .keep_all = TRUE)
    
    plotHistogram(virome = filtered, variable = histSelectedVariable, title = '')  # Use the selected variable for plotting
  })
  
  # TODO: Make this more useful
  # ======================= Scatter Plot ============================
  # output$ScatterPlot <- renderPlotly({
  #   req(plotData())  # Ensure data is available
  #   
  #   scatterXVariable <- input$scatterXVariableSelect 
  #   scatterYVariable <- input$scatterYVariableSelect
  #   req(scatterXVariable) 
  #   req(scatterYVariable)  
  #   
  #   plotScatter(virome = plotData(), x = scatterXVariable, y = scatterYVariable, title = '')  # Use the selected variables for plotting
  # })
  
  # ======================= Violin Plot ============================
  output$ViolinPlot <- renderPlotly({
    req(plotData())  # Ensure data is available
    
    violinXVariable <- input$violinXVariableSelect  
    violinYVariable <- input$violinYVariableSelect 
    
    req(violinXVariable)
    req(violinYVariable)
    
    plotViolin(virome = plotData(), x = violinXVariable, y = violinYVariable, title = '')
  })
  
  # ======================= Virus Positive Plot ============================
  output$plotVirusPositive <- renderPlotly({
    req(viromeData())
    req(serratusProcessed())
    plotVirusPositive(viromeData = viromeData(), processedRunData = serratusProcessed())
        
  })
  
  # ======================= Sankey ============================
  output$sankeyPlot <- renderPlotly({
    req(plotData())
    plotSankey(virome = plotData())
  })
  
  # ======================= Network ============================
  output$viromeNetwork <- renderSimpleNetwork({
    req(viromeData())
    
    networkSplitVariable <- input$networkSplitVariable
    networkCutoff <- input$networkCutoff

    plotViromeNetwork(virome = viromeData(), splitVariable = networkSplitVariable, threshold = networkCutoff)
  })
  
  # ======================= Venn Diagram ============================
  output$vennDiagram <- renderPlot({
    req(viromeData())
    req(vennSet1())
    req(vennSet2())
    req(vennSet3())
    req(vennSet4())
    
    plotsOTUVenn(virome = viromeData(), set1 = vennSet1(), set2 = vennSet2(), set3 = vennSet3(), set4 = vennSet4())
  })
  
  # ======================= PCA ============================
  output$plotPCA <- renderPlotly({
    req(viromeData())
    
    pcaVariable <- input$pcaVariable
    pcaType <- input$pcaType
    
    plotViromePCA(virome = viromeData(), mode = pcaVariable, sOTU = as.logical(pcaType))
  })
}



shinyApp(ui = ui, server = server)

# [END]
