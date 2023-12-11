# test_hostViralBurden.R

library(testthat)

con <- palmid::SerratusConnect()

# ========= getHostViralBurden =========
test_that("getHostViralBurden throws an error when virome is not provided", {
  expect_error(getHostViralBurden())
})

test_that("getHostViralBurden returns a data frame with the correct columns", {
  virome <- getVirome(tax = "Meloidogyne", con = con)
  vbDF <- getHostViralBurden(virome = virome)
  expect_true(is.data.frame(vbDF))
  expect_true('scientific_name' %in% names(vbDF))
  expect_true('total' %in% names(vbDF))
  expect_true('virus_positive' %in% names(vbDF))
  expect_true('percent_virus_positive' %in% names(vbDF))
})

test_that("getHostViralBurden calculates the correct values", {
  virome <- getVirome(tax = "Meloidogyne", con = con)
  vbDF <- getHostViralBurden(virome = virome)
  expect_true(all(vbDF$total >= 0))
  expect_true(all(vbDF$virus_positive >= 0))
  expect_true(all(vbDF$percent_virus_positive >= 0 & vbDF$percent_virus_positive <= 1))
})

# [END]