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
```

## References
A. Babaian and R. C. Edgar (2022), Ribovirus classification by a polymerase barcode sequence. *PeerJ*

Andersen, K. S., Kirkegaard, R. H., Karst, S. M. & Albertsen, M. ampvis2: an R package to analyse and visualise 16S rRNA amplicon data. 299537 Preprint at https://doi.org/10.1101/299537 (2018).

Chamberlain, S., Arendsee, Z. & Stirling, T. taxizedb: Tools for Working with ‘Taxonomic’ Databases. (2023).

R. C. Edgar et al. (2021), Petabase-scale sequence alignment catalyses viral discovery. *Nature*.

Wickham et al., (2019). Welcome to the tidyverse. *Journal of Open Source Software*, 4(43), 1686, https://doi.org/10.21105/joss.01686

Xia, Y. & Sun, J. Bioinformatic and Statistical Analysis of Microbiome Data: From Raw Sequences to Advanced Modeling with QIIME 2 and R. *Springer Nature*, 2023.

