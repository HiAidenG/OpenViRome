## OpenViRome
                                           

An R package for viromics-based analysis of Serratus data. 

### Description

OpenViRome is an R package designed to generate 'virome reports' of Serratus
data given a taxonomic group of interest. Its primary focus is on generating
meaningful and interpretable plots of viral diversity for streamlined analysis,
not the generation of new data. It is meant to be used in tandem with palmid,
which provides a contrasting 'virus-first' approach. OpenViRome can take as
input a taxonomic group of interest, or a vector of SRA run ids that have been
queried by Serratus.

The ultimate goal of the two packages is to create a framework for prioritizing
novel viruses relevant to human health and agriculture. 

### Installation

```r
devtools::install_github("HiAidenG/OpenViRome", build_vignettes = TRUE)
require("OpenViRome")
```

### Overview

```r
ls("OpenViRome")
data("TylenchoideaVirome", package = "OpenViRome")
browseVignettes("OpenViRome")
```

### Contributions

Conceptualization and code for figures: Aiden Hiller. Supervision from Dr. Artem Babaian

Taxonomic lookup made possible with [Taxize](https://github.com/ropensci/taxize)

### References

A. Babaian and R. C. Edgar (2022), Ribovirus classification by a polymerase barcode sequence. PeerJ

R. C. Edgar et al. (2021), Petabase-scale sequence alignment catalyses viral discovery. Nature.

Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686

### Acknowledgements

This package was developed as part of an assessment for 2023 BCB410H: Applied Bioinformat-
ics course at the University of Toronto, Toronto, CANADA. OpenViRome welcomes issues,
enhancement requests, and other contributions. To submit an issue, use the GitHub issues.





