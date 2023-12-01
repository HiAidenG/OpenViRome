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
          fluidRow(
            column(6, plotlyOutput("PieChart")),
            column(6, plotlyOutput("VirusBarChart"))
          ),
          fluidRow(
            column(12, plotlyOutput("Sankey"))
          ),
          fluidRow(
            column(6, plotlyOutput("AlphaDiversity"))
            # column(6, plotlyOutput("BetaDiversityOrdination")) - commented out as per your original code
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

    # Render the plot
    output$PieChart <- renderPlotly({
      req(viromeData()) # Ensure virome data is not NULL

      # Generate and render the pie chart
      plotVirome(virome = viromeData())
    })

    output$VirusBarChart <- renderPlotly({
      req(viromeData()) # Ensure virome data is not NULL

      # Generate and render the bar chart
      plotVirusPositive(getHostViralBurden(virome = viromeData()))
    })

    output$Sankey <- renderPlotly({
      req(viromeData()) # Ensure virome data is not NULL

      # Generate and render the sankey diagram
      drawVirusSankey(virome = viromeData())
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
