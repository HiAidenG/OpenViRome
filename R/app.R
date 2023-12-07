library(shiny)
library(palmid)
#library(ampvis2)
library(plotly)
library(tidyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(taxizedb)
require(circlize)
require(ComplexHeatmap)

# Define UI
ui <- fluidPage(
  titlePanel("openViRome"),

  sidebarLayout(
    sidebarPanel(
      textInput("taxonomicName", "Enter Taxonomic Name:", value = ""),
      actionButton("submit", "Submit"),

      # Add a conditional panel to download the virome
      conditionalPanel(
        condition = "output.stats", # i.e. once virome data is available
        downloadButton("downloadVirome", "Download Virome CSV")
      )

    ),

    mainPanel(
      tabsetPanel(
        # First Tab: All Objective Plots
        tabPanel("Objective Plots",

                 # Add some summary statistics at the top
                 verbatimTextOutput("stats"),

                 fluidRow(
                   column(6, plotlyOutput("PieChart")),
                   column(6, plotlyOutput("VirusBarChart"))
                 ),
                 fluidRow(
                   column(12, plotlyOutput("Sankey"))
                 ),
                 fluidRow(
                   column(12, uiOutput("phylumFilterPanel"))
                 ),
                 fluidRow(
                   column(6, plotlyOutput("AlphaDiversity"))
                   # column(6, plotlyOutput("BetaDiversityOrdination"))
                 )
        ),
        tabPanel("Subjective Plots",
                 plotOutput("viroHeatMap"),
                 plotlyOutput("palmPrevalence")
        )
      )
    )
  )
)


# Define server logic
server <- function(input, output, session) {

  viromeData <- reactiveVal(NULL)
  selectedPhylum <- reactiveVal(NULL)
  # ampvisData <- reactiveVal(NULL)

  # Observe the submit button
  observeEvent(input$submit, {
    con <- palmid::SerratusConnect()
    req(input$taxonomicName) # Ensure input is not empty
    virome <- getVirome(tax = input$taxonomicName, con = con)

    viromeData(virome)

    # TODO: Fix beta diversity
    #ampvisData <- virometoOTUTable(virome = virome)
    #ampvisMeta <- virometoMetadataTable(virome = virome)

    # Generate the ampvis object
    # ampvisObj <- ampvis2::amp_load(otutable = ampvisData, metadata = ampvisMeta)

    # Store data

    # ampvisData(ampvisObj)

    # Calculate summary statistics
    output$stats <- renderPrint({
      req(viromeData()) # Ensure virome data is not NULL

      stats <- getViromeSummary(virome = viromeData())
      # Print virome stats
      cat("Number of unique sOTUs: ", stats[[1]], "\n")
      cat("Mean normalized sOTU coverage: ", stats[[2]], "\n")
      cat("Median normalized sOTU coverage: ", stats[[3]], "\n")
      cat("Max normalized sOTU coverage: ", stats[[4]], "\n")
      cat("Virus positive runs: ", stats[[5]], "\n")
      cat("Runs processed: ", stats[[6]], "\n")

    })

    output$downloadVirome <- downloadHandler(
      filename = function() {
        paste0(input$taxonomicName, "_virome.csv")
      },
      content = function(file) {
        write.csv(virome[[1]], file, row.names = FALSE)
      }
    )

    output$PieChart <- renderPlotly({
      req(viromeData()) # Ensure virome data is not NULL
      plotViromePie(virome = viromeData())
    })

    output$Sankey <- renderPlotly({
      req(viromeData()) # Ensure virome data is not NULL
      current_phylum <- selectedPhylum()

      if (is.null(current_phylum) || current_phylum == "All") {
        drawVirusSankey(virome = viromeData())
      }
      else {
        # TODO: Very hacky solution to match colors, but works for now
        phyla <- getAvailablePhyla(virome = viromeData())
        if (length(phyla) > 9) {
          colors <- RColorBrewer::brewer.pal(n = length(phyla),
                                             name = "Set3")
        }
        else {
          colors <- grDevices::palette.colors(n = length(phyla),
                                              palette = "Okabe-Ito",
                                              alpha = 0.7)
        }
        names(colors) <- phyla
        color <- colors[current_phylum]

        drawVirusSankey(virome = viromeData(), phylumFilter = current_phylum,
                        colors = color)
      }
    })

    output$phylumFilterPanel <- renderUI({
      req(viromeData()) # Only render this UI if viromeData is available

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



    output$VirusBarChart <- renderPlotly({
      req(viromeData()) # Ensure virome data is not NULL

      # Generate and render the bar chart
      plotVirusPositive(getHostViralBurden(virome = viromeData()))
    })


    output$AlphaDiversity <- renderPlotly({
      req(viromeData()) # Ensure virome data is not NULL

      # Generate and render the alpha diversity plot
      plotAlphaDiversity(virome = viromeData())
    })

    # output$BetaDiversityOrdination <- renderPlotly({
      # req(viromeData()) # Ensure virome data is not NULL

      # Generate and render the beta diversity plot
      # ampvis2::amp_ordinate(ampvisData(), sample_color_by="scientific_name",
      #                       filter_species = 0.0,
      #                       sample_plotly = "all")
    # })
  })

  output$viroHeatMap <- renderPlot({
    req(viromeData()) # Ensure virome data is not NULL

    # Generate and render the heatmap
    viroMap(virome = viromeData(), minCov = 1)
  })

  output$palmPrevalence <- renderPlotly({
    req(viromeData()) # Ensure virome data is not NULL

    # Generate and render the prevalence plot
    palmPrevalence(virome = viromeData())
  })
}

# Run the application
shinyApp(ui = ui, server = server)
