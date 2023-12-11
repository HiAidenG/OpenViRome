library(testthat)
library(dplyr)

# Load test data
data("TylenchoideaVirome")

# Test getAlphaDiversity with default parameters
test_that("getAlphaDiversity calculates Shannon diversity by default", {
  result <- getAlphaDiversity(virome = TylenchoideaVirome)
  expect_true("shannon" %in% colnames(result))
  expect_false("simpson" %in% colnames(result))
  expect_false("richness" %in% colnames(result))
  expect_false("evenness" %in% colnames(result))
})

# Test getAlphaDiversity with specific bioSample
test_that("getAlphaDiversity calculates Shannon diversity for a specific bioSample", {
  result <- getAlphaDiversity(virome = TylenchoideaVirome, bioSample = "SAMN25258043")
  expect_true("shannon" %in% colnames(result))
  expect_false("simpson" %in% colnames(result))
  expect_false("richness" %in% colnames(result))
  expect_false("evenness" %in% colnames(result))

  # should only have the one biosample
  expect_equal(nrow(result), 1)
  expect_equal(result$bio_sample[1], "SAMN25258043")
})

# Test getAlphaDiversity with mode = "simpson"
test_that("getAlphaDiversity calculates Simpson diversity", {
  result <- getAlphaDiversity(virome = TylenchoideaVirome, mode = "simpson")
  expect_false("shannon" %in% colnames(result))
  expect_true("simpson" %in% colnames(result))
  expect_false("richness" %in% colnames(result))
  expect_false("evenness" %in% colnames(result))
})

# Test getAlphaDiversity with mode = c("shannon", "simpson", "richness", "evenness")
test_that("getAlphaDiversity calculates all metrics", {
  result <- getAlphaDiversity(virome = TylenchoideaVirome, mode = c("shannon", "simpson", "richness", "evenness"))
  expect_true("shannon" %in% colnames(result))
  expect_true("simpson" %in% colnames(result))
  expect_true("richness" %in% colnames(result))
  expect_true("evenness" %in% colnames(result))
})

# Test getAlphaDiversity with invalid mode
test_that("getAlphaDiversity throws an error for invalid mode", {
  expect_error(getAlphaDiversity(virome = TylenchoideaVirome, mode = "invalid_mode"))
})

# Test getAlphaDiversity without providing a virome object
test_that("getAlphaDiversity throws an error without a virome object", {
  expect_error(getAlphaDiversity())
})

# Test getAlphaDiversity calculates metrics correctly
test_that("getAlphaDiversity calculates alpha metrics correctly" , {
  result <- getAlphaDiversity(virome = TylenchoideaVirome, mode = c("shannon", "simpson", "richness"))
  expect_equal(result$shannon[5], 0)
  expect_equal(result$simpson[5], 1)
  expect_equal(result$richness[5], 1)
})
