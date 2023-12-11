# app.R

library(shiny)
library(palmid)
library(plotly)
library(tidyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(bs4Dash)
require(circlize)
require(grDevices)
require(fresh)
require(OpenViRome)

viromeTheme <- fresh::create_theme(
  bs4dash_status(
    primary = "#6495ED",
    secondary = "#818589",
    success = "#00CD66",
  )
)

ui <- bs4DashPage(
  freshTheme = viromeTheme,
  title = "openViRome",
  header = bs4Dash::dashboardHeader(
    title = "openViRome v0.02",
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
        text = "Virome Diversity",
        tabName = "ViromeDiversity",
        icon = icon("chart-bar")
      ),
      # menuItem(
      #   text = "Source-sOTU Heatmap",
      #   tabName = "sOTUHeatmap",
      #   icon = icon("vial-virus")
      # ),
      menuItem(
        text = "SRA Prevalence",
        tabName = "sOTUScatter",
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
            h2("openViRome"),
            p("Explore Serratus data through a viromics lens! openViRome aims to provide a user-friendly interface for virus discovery by democratizing public sequencing data.")
          )
        ), # Row 1 close
        fluidRow(
          bs4Dash::box(title = "Query Serratus",
                       p("Here you can query Serratus data either by taxonomic name (at any rank, like 'Meloidogyne') or by comma-separated SRA accessions. Be aware that large taxonomic groups may take a long time to query."),
                       closable = FALSE,
                       collapsible = FALSE,
                       width = 6,
                       solidHeader = TRUE,
                       status = "primary",
                       textInput("taxonomicName", label = NULL,
                                 placeholder = "Taxonomic Name",
                                 width = '50%'),

                       p("OR"),
                       textInput("accessions", label = NULL,
                                 placeholder = "SRA Accessions (comma separated)",
                                 width = '50%'),
                       actionButton("submitQuery", "Submit", style = "width: 50%;")
                       ),
          bs4Dash::box(title = "Upload Virome",
                                p("Upload a formatted virome file in CSV format."),
                                p("See the 'help' tab for more information."),
                                p("Example .csv provided at OpenViRome/inst/extdata/Meloidogyne_virome.csv"),
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
          p("This tool is currently still in beta. Please report any bugs or issues to the openVirome GitHub page."),
          h4("How to use openViRome"),
          p("openViRome is designed to be able to provide a 'virome-style' analysis for any arbitrary set of SRA data. As of v0.02, openViRome can be used in three ways:"),
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
              width = 6,
              uiOutput("VirusBarChartCard")
            ),
            column(
              width = 6,
              uiOutput("PieChart")
            )
        ),
        fluidRow(
          column(
            width = 12,
            uiOutput("SankeyCard")
          )
        ),
        fluidRow(
          column(
            width = 6,
            uiOutput("phylumFilterPanel")
          )
        )
      ), # Close ViromeOverview
      bs4TabItem(
        tabName = "ViromeDiversity",
        h2("Virome Diversity"),
        fluidRow(
          column(
            width = 12,
            uiOutput("AlphaDiversityCard")
          )
        ),
        fluidRow(
          bs4Dash::bs4Card(
            title = "Beta Diversity",
            plotlyOutput("BetaDiversityPlot", height = "600px"),
            width = 12,
            collapsible = FALSE
          )
        )
      ), # Close ViromeDiversity
      # TODO: heatmap still doesn't display properly
      # bs4TabItem(
      #   tabName = "sOTUHeatmap",
      #   h2("sOTU Heatmap"),
      #   fluidRow(
      #     bs4Dash::bs4Card(
      #       title = "Heatmap",
      #       plotOutput("viroHeatMap", height = "600px"),
      #       width = 12,
      #       collapsible = TRUE,
      #       maximizable = TRUE
      #     )
      #   )
      # ),
      bs4TabItem(
        tabName = "sOTUScatter",
        h2("sOTU Scatterplot"),
        fluidRow(
          bs4Dash::bs4Card(
            title = "sOTU Scatterplot",
            plotlyOutput("palmPrevalence", height = "600px"),
            width = 12,
            collapsible = FALSE
          )
        )
        # TODO: Add distributions for the three variables
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
  selectedPhylum <- reactiveVal(NULL)

  # ======================= Query button ============================
  observeEvent(input$submitQuery, {
    showModal(modalDialog(
      title = "Processing...",
      "Your data is being loaded. Depending on the size of your query, this may take a few minutes.",
      easyClose = TRUE,
      footer = NULL
    ))

    con <- palmid::SerratusConnect()
    errorOccurred <- FALSE

    if (!is.null(input$taxonomicName)) {
      req(input$taxonomicName) # Ensure input is not empty
      tryCatch({
        virome <- getVirome(tax = input$taxonomicName, con = con)
      }, error = function(e) {
        showNotification("Error: There was an issue loading your virome. Check your taxon for typos.", type = "error")
        errorOccurred <<- TRUE
      }, finally = {
        removeModal()
      })
    }
    else if (!is.null(input$accessions)) {
      req(input$accessions) # Ensure input is not empty
      tryCatch({
        virome <- getVirome(sra = input$accessions, con = con)
      }, error = function(e) {
        showNotification("Error: There was an issue loading your virome. Are you sure your accessions are comma-separated?", type = "error")
        errorOccurred <<- TRUE
      }, finally = {
        removeModal()
      })
    }

    if (!errorOccurred) {
      viromeData(virome)
      removeModal()
    }
  })


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

      con <- palmid::SerratusConnect()

      # TODO: this doesn't work, but it's a temporary fix
      # Use first row as tax
      tax <- virome[1, "scientific_name"][[1]]
      tax <- strsplit(tax, " ")[[1]][1]
      runDF <- taxLookup(tax = tax, con = con)
      runDF <- runDF %>%
        dplyr::mutate(virus_positive = ifelse(run %in% virome$run, TRUE, FALSE))
      virome <- list(virome, runDF)

      viromeData(virome)

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
      filename = function() {
        if (is.null(input$submitQuery)) {
          return(paste0("virome_", input$uploadFile$name))
        }
        else {
          if (!is.null(input$taxonomicName)) {
            return(paste0(input$taxonomicName, "_virome.csv"))
          }
          else if (!is.null(input$accessions)) {
            return(paste0(input$accessions, "_virome.csv"))
          }
        }
      },
      content = function(file) {
        data <- viromeData()
        toWrite <- data[[1]]
        write.csv(toWrite, file = file, row.names = FALSE)
      }
    )

    # ======================= Summary stats ============================
    output$statsCard <- renderUI({

      if (!is.null(viromeData())) {
      stats <- getViromeSummary(virome = viromeData())

      # Create a bs4Card to display the stats
      bs4Card(
        title = "Virome Summary Stats",
        status = "success",
        solidHeader = TRUE,
        collapsible = FALSE,
        width = NULL,
        p("Number of unique sOTUs: ", stats[[1]]),
        p("Mean normalized sOTU coverage: ", round(stats[[2]], 2)),
        p("Median normalized sOTU coverage: ", round(stats[[3]], 2)),
        p("Max normalized sOTU coverage: ", round(stats[[4]], 2)),
        p("Virus positive runs: ", stats[[5]]),
        p("Runs processed: ", stats[[6]])
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

    # ======================= Pie Chart ============================
    output$PieChart <- renderUI({
      req(viromeData()) # Ensure data is available
      bs4Dash::box(
        title = "Pie Chart",
        status = "primary",
        solidHeader = FALSE,
        width = NULL,
        collapsible = FALSE,
        plotViromePie(virome = viromeData())
      )
    })

    # ======================= Sankey ================================
    output$Sankey <- renderPlotly({
      req(viromeData())
      current_phylum <- selectedPhylum()

      if (is.null(current_phylum) || current_phylum == "All") {
        p <- drawVirusSankey(virome = viromeData())
      } else {
        # very hacky solution to match colors to phyla, but works for now.
        phyla <- getAvailablePhyla(virome = viromeData())
        if (length(phyla) > 9) {
          colors <- RColorBrewer::brewer.pal(n = length(phyla), name = "Set3")
        } else {
          colors <- grDevices::palette.colors(n = length(phyla),
                                              palette = "Okabe-Ito",
                                              alpha = 0.7)
        }
        names(colors) <- phyla
        color <- colors[current_phylum]

        p <- drawVirusSankey(virome = viromeData(), phylumFilter = current_phylum,
                        colors = color)
      }
      p <- p %>% layout(autosize = TRUE)

      return(p)
    })

    output$SankeyCard <- renderUI({
      bs4Dash::box(
        title = "Virome Sankey Diagram",
        status = "primary",
        solidHeader = FALSE,
        collapsible = FALSE,
        maximizable = TRUE,
        width = NULL,
        plotlyOutput("Sankey", width = "100%", height = "100%")
      )
    })

    output$phylumFilterPanel <- renderUI({
      req(viromeData())

      availablePhyla <- getAvailablePhyla(virome = viromeData())

      # Create a selectInput for phylum filtering
      selectInput("selectedPhylum", "Select Phylum:",
                  choices = c("All", availablePhyla),
                  selected = "All")
    })

    # Update the selectedPhylum reactive value whenever the input changes
    observe({
      selectedPhylum(input$selectedPhylum)
    })

    # ======================= Bar Chart ================================
    output$VirusBarChart <- renderPlotly({
      req(viromeData())

      # Generate and render the bar chart
      plotVirusPositive(getHostViralBurden(virome = viromeData()))
    })

    output$VirusBarChartCard <- renderUI({
      bs4Dash::box(
        title = "Virus Positive Runs",
        status = "primary",
        solidHeader = FALSE,
        collapsible = FALSE,
        maximizable = FALSE,
        width = NULL,
        plotlyOutput("VirusBarChart", width = "100%", height = "100%")
      )
    })

    # ======================= Alpha Diversity ============================
    output$AlphaDiversity <- renderPlotly({
      req(viromeData())

      plotAlphaDiversity(virome = viromeData())
    })

    output$AlphaDiversityCard <- renderUI({
      bs4Dash::box(
        title = "BioSample Alpha Diversity",
        status = "primary",
        solidHeader = FALSE,
        collapsible = FALSE,
        maximizable = FALSE,
        width = NULL,
        plotlyOutput("AlphaDiversity", width = "100%", height = "100%")
      )
    })

    # ======================= Beta Diversity ============================
    output$BetaDiversityPlot <- renderPlotly({
        req(viromeData())
        plotBetaDiversity(virome = viromeData())
      })

  # ======================= Heatmap ====================================
    output$viroHeatMap <- renderPlot({
      req(viromeData())
      heatmapPlot <- viroMap(virome = viromeData(), minCov = 1)
    })

  # ======================= Scatter plot ============================
    output$palmPrevalence <- renderPlotly({
      req(viromeData())
      palmPrevalence(virome = viromeData())
    })

    # output$sradist <- renderPlotly({
    #   req(viromeData())
    #
    # # plot the distribution of sra runs
    #   virome <- viromeData()
    #   data <- virome[[1]]
    #   sotus <- data %>% select(sotu) %>% unique()
    # })

  # TODO: Add plots for 4 distributions: sra runs, bioprojects, coverage, and genbank identity
}

shinyApp(ui = ui, server = server)

# [END]
