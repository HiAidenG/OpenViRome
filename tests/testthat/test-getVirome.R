# Test for getVirome

library(testthat)
library(OpenViRome)

test_that("getVirome throws an error when both tax and sra are NULL", {
  con <- palmid::SerratusConnect()
  expect_error(getVirome(con = con))
})

test_that("getVirome throws an error when both tax and sra are provided", {
  con <- palmid::SerratusConnect()
  expect_error(getVirome(tax = "TaxTerm", sra = c("SRR12345"), con = con))
})

test_that("getVirome throws an error when no connection is provided", {
  expect_error(getVirome(tax = "TaxTerm"))
})



# [END]
