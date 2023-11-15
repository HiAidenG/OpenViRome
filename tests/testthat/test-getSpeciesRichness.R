# Tests for getSpeciesRichness function

library(testthat)
require(OpenViRome)

test_that("getSpeciesRichness returns correct value", {
  virome <- data.frame(scientific_name = c("Species1", "Species2"),
                       tax_family = c("Family1", "Family2"))
  expect_equal(getSpeciesRichness(virome), 2)
})

test_that("getSpeciesRichness handles sourceSpecies filtering", {
  virome <- data.frame(scientific_name = c("Species1", "Species2", "Species1"),
                       tax_family = c("Family1", "Family2", "Family1"))
  expect_equal(getSpeciesRichness(virome, "Species1"), 1)
})

test_that("getSpeciesRichness throws error for invalid sourceSpecies", {
  virome <- data.frame(scientific_name = c("Species1", "Species2"),
                       tax_family = c("Family1", "Family2"))
  expect_error(getSpeciesRichness(virome, "Species3"))
})

# [END]
