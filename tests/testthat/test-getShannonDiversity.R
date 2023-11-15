# Tests for getShannonDiversity function

library(testthat)
require(OpenViRome)

test_that("getShannonDiversity returns correct value", {
  virome <- data.frame(scientific_name = c("Species1", "Species2"),
                       tax_family = c("Family1", "Family2"))
  expect_equal(getShannonDiversity(virome), -(log(1/2) * 1/2 * 2))
})

test_that("getShannonDiversity handles sourceSpecies filtering", {
  virome <- data.frame(scientific_name = c("Species1", "Species2", "Species1"),
                       tax_family = c("Family1", "Family2", "Family1"))
  expect_equal(getShannonDiversity(virome, "Species1"), 0)
})

test_that("getShannonDiversity throws error for invalid sourceSpecies", {
  virome <- data.frame(scientific_name = c("Species1", "Species2"),
                       tax_family = c("Family1", "Family2"))
  expect_error(getShannonDiversity(virome, "Species3"))
})

# [END]
