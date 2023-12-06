library(shiny)
library(palmid)
library(ampvis2)
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
      actionButton("submit", "Submit")
    ),

    mainPanel(
      # Using tabsetPanel to create tabs
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
            column(6, plotlyOutput("AlphaDiversity"))
            # column(6, plotlyOutput("BetaDiversityOrdination"))
          )
        ),
        # Second Tab: Blank for now
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

  # Reactive value for virome data
  viromeData <- reactiveVal(NULL)
  selectedPhylum <- reactiveVal(NULL)
  # ampvisData <- reactiveVal(NULL)

  # Observe the submit button
  observeEvent(input$submit, {
    con <- palmid::SerratusConnect()
    req(input$taxonomicName) # Ensure input is not empty
    virome <- getVirome(tax = input$taxonomicName, con = con)

    # TODO: Fix beta diversity
    #ampvisData <- virometoOTUTable(virome = virome)
    #ampvisMeta <- virometoMetadataTable(virome = virome)

    # Generate the ampvis object
    # ampvisObj <- ampvis2::amp_load(otutable = ampvisData, metadata = ampvisMeta)

    # Store data
    viromeData(virome)
    # ampvisData(ampvisObj)

    # Calculate summary statistics
    output$stats <- renderPrint({
      req(viromeData()) # Ensure virome data is not NULL

      stats <- getViromeSummary(virome = viromeData())
      numSOTUs <- stats['numSOTUs'][[1]]
      medCoverage <- stats['medianNormCov'][[1]]
      maxCoverage <- stats['maxNormCov'][[1]]

      # Return the summary statistics
      paste("Number of Unique sOTUs: ", numSOTUs, "\n",
            "Median Normalized Coverage ", medCoverage, "\n",
            "Max Normalized Coverage: ", maxCoverage, "\n",
            sep = "")
    })


    output$PieChart <- renderPlotly({
    req(viromeData()) # Ensure virome data is not NULL

    pie_data <- viromeData() %>%
                distinct(sotu) %>%
                group_by(tax_phylum) %>%
                summarise(Count = n()) %>%
                arrange(desc(Count))

    # Create the pie chart directly here
    p <- plot_ly(pie_data, labels = ~tax_phylum, values = ~Count, type = 'pie', 
                 key = ~Phylum, source = "PieChart") %>%
         layout(title = 'Distribution by Phylum')

    # Register the 'plotly_click' event
    event_register(p, 'plotly_click')

    p # Return the plot
})

observeEvent(event_data("plotly_click", source = "PieChart"), function(click_data) {
    selectedPhylum(click_data$key)
})


    output$Sankey <- renderPlotly({
        req(viromeData()) # Ensure virome data is not NULL
        current_phylum <- selectedPhylum()

        drawVirusSankey(virome = viromeData(), phylumFilter = current_phylum)
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
