# Tests for getSimpsonDiversity function

library(testthat)
require(OpenViRome)

test_that("getSimpsonDiversity returns correct value", {
  virome <- data.frame(scientific_name = c("Species1", "Species2"),
                       tax_family = c("Family1", "Family2"))
  expect_equal(getSimpsonDiversity(virome), 1/2)
})

test_that("getSimpsonDiversity handles sourceSpecies filtering", {
  virome <- data.frame(scientific_name = c("Species1", "Species2", "Species1"),
                       tax_family = c("Family1", "Family2", "Family1"))
  expect_equal(getSimpsonDiversity(virome, "Species1"), 1)
})

test_that("getSimpsonDiversity throws error for invalid sourceSpecies", {
  virome <- data.frame(scientific_name = c("Species1", "Species2"),
                       tax_family = c("Family1", "Family2"))
  expect_error(getSimpsonDiversity(virome, "Species3"))
})

# [END]
