# OpenViRome
An R package for querying and interpreting Serratus data through a viromics lens.

## Description:
OpenViRome is an R package designed to generate 'virome reports' of Serratus
data given a taxonomic group of interest, or any arbitrary set of SRA runs.
Its primary focus is on generating meaningful and interpretable plots of viral 
diversity for streamlined analysis, not the generation of new data. It is meant
to be used in tandem with palmid,which provides a contrasting 'virus-first' 
approach. OpenViRome can take as input a taxonomic group of interest, or a 
vector of SRA run ids provided they have been processed by Serratus. 

OpenViRome provides functionality for querying the Serratus database, understanding
the taxonomic composition of the virome, calculating alpha and beta diversity
metrics, and visualizing the results. It also provides functions for plotting
the prevalence of viral sOTUs across the SRA. 

Most of the functionality of OpenViRome is available through the Shiny app,
which generates no-code summary reports for a user provided input. 

## Installation:
To install the latest release:
```{r}
install.packages("devtools")
library("devtools")
devtools::install_github("HiAidenG/OpenViRome", build_vignettes = TRUE)
require("OpenViRome")
```

To run the shiny app:
```{r}
OpenViRome::viromeDashboard()
```

## Overview:
```{r}
ls("package:OpenViRome")
data(package = "OpenViRome")
browseVignettes("OpenViRome")
```

The main function of OpenViRome is `getVirome()`, which takes as input either a
taxonomic term or a vector of SRA run ids. It returns a 'virome' object, which
is a list of two dataframes containing the virome data and the run data.

A few of the functions in OpenViRome require a connection to the Serratus
database. This can be established using the `SerratusConnect()` function from
the palmid package.

Some of the other functions include: 
- `plotViromePie()` - plots the taxonomic composition of the virome as a pie chart
- `drawVirusSankey()` - plots the taxonomic composition of the virome and the source bioproject of    each viral sOTU as a sankey diagram.
- `plotBetaDiversity()` - plots bray-curtis dissimilarity using PCoA via ampvis2. Functions for      converting a virome object to an ampvis2 object are also provided.
- `plotAlphaDiversity()` - plots alpha diversity metrics (e.g. Shannon, Simpson).

Additional information and functions are covered in the vignette.

![](./inst/extdata/openViRome.png)

## Contributions:
All functions were authored by Aiden Hiller. Some of the functions depend on packages written by others, including:
 - `ampvis2` - for beta diversity analysis (`plotBetaDiversity()`).
 - `palmid` - for querying the Serratus database. 
 - `taxizedb` - for querying NCBI taxonomy (`taxLookup()`).
 
Additionally, `ggplot` and `plotly` are used for plotting throughout the package.

The shiny app was written by Aiden Hiller, with help from the `shinydashboard` and `bs4dash` packages.

Provided data (`TylenchoideaVirome`) was generated as part of the Serratus project (Edgar et al., 2022). 

## References
A. Babaian and R. C. Edgar (2022), Ribovirus classification by a polymerase barcode sequence. *PeerJ*

Andersen, K. S., Kirkegaard, R. H., Karst, S. M. & Albertsen, M. ampvis2: an R package to analyse and visualise 16S rRNA amplicon data. 299537 Preprint at https://doi.org/10.1101/299537 (2018).

Chamberlain, S., Arendsee, Z. & Stirling, T. taxizedb: Tools for Working with ‘Taxonomic’ Databases. (2023).

R. C. Edgar et al. (2021), Petabase-scale sequence alignment catalyses viral discovery. *Nature*.

Wickham et al., (2019). Welcome to the tidyverse. *Journal of Open Source Software*, 4(43), 1686, https://doi.org/10.21105/joss.01686

Xia, Y. & Sun, J. Bioinformatic and Statistical Analysis of Microbiome Data: From Raw Sequences to Advanced Modeling with QIIME 2 and R. *Springer Nature*, 2023.

## Acknowledgements

This package was developed as part of an assessment for 2023 BCB410H: Applied Bioinformat-
ics course at the University of Toronto, Toronto, CANADA. OpenViRome welcomes issues,
enhancement requests, and other contributions. To submit an issue, use the GitHub issues.
