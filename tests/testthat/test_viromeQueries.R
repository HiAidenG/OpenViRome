# test_viromeQueries.R

library(testthat)

con <- palmid::SerratusConnect()

# ========= getVirome =========
test_that("getVirome returns a virome object when tax is provided", {
  virome <- getVirome(tax = "Meloidogyne", con = con)
  expect_true(is.list(virome))
  expect_true(length(virome) == 2)
  expect_true(is.data.frame(virome[[1]]))
  expect_true('run' %in% names(virome[[1]]))
  expect_true(is.data.frame(virome[[2]]))
})

test_that("getVirome returns a virome object when sra is provided", {
  virome <- getVirome(sra = c("SRR17756040", "SRR5942326"), con = con)
  expect_true(is.list(virome))
  expect_true(length(virome) == 2)
  expect_true(is.data.frame(virome[[1]]))
  expect_true('run' %in% names(virome[[1]]))
  expect_true(is.data.frame(virome[[2]]))
})

test_that("getVirome throws an error when con is not provided", {
  expect_error(getVirome(tax = "Meloidogyne"))
})

test_that("getVirome throws an error when neither tax nor sra is provided", {
  expect_error(getVirome(con = con))
})

test_that("getVirome throws an error when both tax and sra are provided", {
  expect_error(getVirome(tax = "Meloidogyne", sra = c("SRR17756040", "SRR5942326"), con = con))
})

test_that("getVirome throws an error when taxonomic term is not found in NCBI taxonomy database", {
  expect_error(getVirome(tax = "NonexistentTerm", con = con))
})

# ========= taxLookup =========
test_that("taxLookup throws an error when con is not provided", {
  expect_error(taxLookup(tax = "Meloidogyne"))
})

test_that("taxLookup throws an error when tax is not provided", {
  expect_error(taxLookup(con = con))
})

test_that("taxLookup throws an error when taxonomic term is not found in NCBI taxonomy database", {
  expect_error(taxLookup(tax = "NonexistentTerm", con = con))
})

test_that("taxLookup returns a data frame of SRA accessions when valid tax and con are provided", {
  result <- taxLookup(tax = "Meloidogyne", con = con)
  expect_true(is.data.frame(result))
  expect_true('run' %in% names(result))
})

# ========= getVirusTaxonomy =========
test_that("getVirusTaxonomy throws an error when otu is not provided", {
  expect_error(getVirusTaxonomy(con = con))
})

test_that("getVirusTaxonomy throws an error when con is not provided", {
  expect_error(getVirusTaxonomy(otu = "VirusName"))
})

test_that("getVirusTaxonomy returns a character vector of taxons when valid otu and con are provided", {
  result <- getVirusTaxonomy(otu = "Flammulina betapartitivirus", con = con)
  expect_true(is.data.frame(result))
  expect_true('tax_phylum' %in% names(result))
})

test_that("getVirusTaxonomy returns the correct taxonomic information when valid otu and con are provided", {
  result <- getVirusTaxonomy(otu = "u16", con = con)
  expect_true(result$tax_phylum == 'Negarnaviricota')
})

# [END]
