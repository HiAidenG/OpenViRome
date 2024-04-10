# basicPlots.R

# This module provides several general functions for plotting data associated with a virome.

#' @title getDefaultViromePalette
#' 
#' @description Returns a default color palette for virome plots (Okabe-ito).
#' 
#' @param numColors The number of colors to return. Default is 9.
#' 
#' @return A vector of colors.
#' 
#' @examples
#' getDefaultViromePalette()
#' getDefaultViromePalette(numColors = 5)
#' 
#' @export
getDefaultViromePalette <- function(numColors = 9) {
    # By default, the okabe ito palette is used.
    colPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")

    if (numColors > length(colPalette)) {
        warning("Number of colors requested exceeds the number of colors in the default palette. Using the default palette.")
        return(colPalette)
    }

    return(colPalette[1:numColors])
}


getViromeVariables <- function(virome = NULL) {
    if (is.null(virome)) {
        stop("No virome data provided.")
    }
    else {
        columns <- colnames(virome)
        # Return all columns in supported params
        categoricalParams <- retrieveParamInfo(paramType = "categoricalPlotting")
        numericalParams <- retrieveParamInfo(paramType = "numericalPlotting")
        supportedParams <- c(categoricalParams, numericalParams)
        validColNames <- intersect(columns, supportedParams)
        return(validColNames)
    }
}

# Helper for sanity checking user-provided arguments to plotting functions
# Checks if the virome object is a tibble, has the specified variables, and the variables are of the correct type for the given plot
checkUserArgs <- function(virome = NULL, variables = NULL, plotType = NULL) {
    if (!tibble::is_tibble(virome)) {
        stop("The virome object must be a tibble.")
    }

    if (is.null(variables)) {
        stop("Please provide a list of variables to plot.")
    }

    if (length(variables) == 0) {
        stop("Please provide at least one variable to plot.")
    }

    if (is.null(plotType)) {
        stop("Please provide a plot type.")
    }

    if (!plotType %in% c("scatter", "histogram", "pie", "violin", "bar")) {
        stop("Invalid plot type. Please provide a valid plot type, one of: 'scatter', 'histogram', 'pie', 'violin', 'bar'.")
    }

    # Check if the variables are in the virome
    viromeVariables <- getViromeVariables(virome)
    if (!all(variables %in% viromeVariables)) {
        stop("One or more variables are not present in the virome.")
    }

    numericVariables <- retrieveParamInfo(paramType = "numericalPlotting")
    categoricalVariables <- retrieveParamInfo(paramType = "categoricalPlotting")
    
    # Check if the variables are of the correct type for the plot
    if (plotType == "scatter") {
        # Both should be numeric
        if (!all(variables %in% numericVariables)) {
            stop("One or more variables are not numeric and cannot be plotted as a scatter plot.")
        }
    } else if (plotType == "histogram") {
        # Only one variable is needed, and it should be numeric
        if (!all(variables %in% numericVariables)) {
            stop("One or more variables are not numeric and cannot be plotted as a histogram.")
        }
    } else if (plotType == "pie") {
        # Only one variable is needed, and it should be categorical
        if (!all(variables %in% categoricalVariables)) {
            stop("One or more variables are not categorical and cannot be plotted as a pie chart.")
        }
    } else if (plotType == "violin") {
        # Two variables are needed, one categorical and one numeric
        if (!all(variables %in% c(numericVariables, categoricalVariables))) {
            stop("One or more variables are not of the correct type for a violin plot.")
        }
    } else if (plotType == "bar") {
        # One variable needed, and it should be categorical
        if (!all(variables %in% categoricalVariables)) {
            stop("One or more variables are not categorical and cannot be plotted as a bar plot.")
        }
    }
}

plotScatter <- function(virome = NULL, x = NULL, y = NULL, title = "Virome Scatter Plot", color = getDefaultViromePalette(numColors = 1), alpha = 0.5, size = 1, jitter = 0.0) {
  checkUserArgs(virome = virome, variables = c(xAxis, yAxis), plotType = "scatter")

  # Create the ggplot scatter plot
  p <- ggplot(virome, aes_string(x = xAxis, y = yAxis)) +
    geom_point(color = color[1], alpha = alpha, size = size) +
    labs(title = title, x = x, y = y) +
    theme_cowplot()

  # Add jitter if specified
  if (jitter > 0) {
    p <- p + geom_jitter(width = jitter)
  }

  plotly_object <- ggplotly(p) %>%
    layout(hovermode = 'closest') %>%
    add_trace(data = virome, x = ~get(xAxis), y = ~get(yAxis),
              type = 'scatter', mode = 'markers',
              marker = list(color = color[1], size = size, opacity = alpha),
              text = ~paste(sOTU, paste(xAxis, ':', get(xAxis)), paste(yAxis, ':', get(yAxis)), sep = '<br>'),
              hoverinfo = 'text')

  return(plotly_object)
}

plotHistogram <- function(virome = NULL, variable = NULL, title = "Histogram with Density and Rug", bins = 30, fill = "#009E73") {
  # Check user arguments
  checkUserArgs(virome = virome, variables = c(variable), plotType = "histogram")

  # Create the histogram with ggplot
  p <- ggplot(virome, aes_string(x = variable)) +
    geom_histogram(aes(y = after_stat(count)), bins = bins, fill = fill) +
    geom_density(alpha = 0.5, fill = fill) +
    geom_rug() +
    labs(title = title, x = variable, y = "Count") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
  plotly_object <- ggplotly(p, tooltip = c("x", "y"))

  return(plotly_object)
}

plotPie <- function(virome = NULL, variable = NULL, title = "Pie Chart", colors = getDefaultViromePalette(numColors = 9)) {
  # Check user arguments
  checkUserArgs(virome = virome, variables = c(variable), plotType = "pie")

  # Determine the number of unique categories
  numCategories <- length(unique(virome[[variable]]))

  # Check if the number of unique categories exceeds the number of colors
  if (numCategories > length(colors)) {
    stop("Number of unique categories exceeds the number of colors in the palette. It's likely that this pie chart will not be visually informative. Try filtering the data first.")
  }

  # Plot directly with plotly, since the ggplot pie charts are not as visually appealing
  plotly_object <- plot_ly(virome, labels = ~get(variable), type = 'pie', marker = list(colors = colors),
                           textinfo = 'label+percent', textposition = 'inside') %>%
    layout(title = title)

  return(plotly_object)
}

plotViolin <- function(virome = NULL, x = NULL, y = NULL, title = "Violin Plot", colors = getDefaultViromePalette(numColors = 9), alpha = 0.5) {
  # Check user arguments
  checkUserArgs(virome = virome, variables = c(x, y), plotType = "violin")
  
  # Figure out which variable is categorical
  if (x %in% retrieveParamInfo(paramType = "categoricalPlotting")) {
    categoricalVariable <- x
    numericVariable <- y
  } else {
    categoricalVariable <- y
    numericVariable <- x
  }
  
  numCategories <- length(unique(virome[[categoricalVariable]]))
  if (numCategories > 9) {
    stop("Number of unique categories exceeds 9. Violin plots are not ideal for visualizing this many categories. Try filtering the data first.")
  }
  
  # Create the violin plot with ggplot
  p <- ggplot(virome, aes_string(x = categoricalVariable, y = numericVariable, fill = categoricalVariable, text = 'sOTU')) +
    geom_violin(alpha = alpha, drop = TRUE) +
    scale_fill_manual(values = colors) +
    geom_jitter(width = 0.1, alpha = alpha, height=0.2) +
    labs(title = title, x = categoricalVariable, y = numericVariable) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill = FALSE)  # Remove the legend
  
  # Convert the ggplot object to a plotly object and specify the tooltip
  plotly_object <- ggplotly(p, tooltip = "text")
  
  return(plotly_object)
}

plotBar <- function(virome = NULL, variable = NULL, title = "Bar Plot", colors = getDefaultViromePalette(numColors = 9)) {
  # Check user arguments
  checkUserArgs(virome = virome, variables = c(variable), plotType = "bar")

  # Determine the number of unique categories
  numCategories <- length(unique(virome[[variable]]))

  # Check if the number of unique categories exceeds the number of colors
  if (numCategories > length(colors)) {
    stop("Number of unique categories exceeds the number of colors in the palette. It's likely that this bar chart will not be visually informative. Try filtering the data first.")
  }

  # Create the bar plot with ggplot
  p <- ggplot(virome, aes_string(x = variable, fill = variable)) +
    geom_bar() +
    scale_fill_manual(values = colors) +
    labs(title = title, x = variable, y = "Count") +
    theme_cowplot() +
    guides(fill = FALSE)

  plotly_object <- ggplotly(p)

  return(plotly_object)
}

# [END]