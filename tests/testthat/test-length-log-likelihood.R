# Test suite for get_length_log_likelihood function
# Tests the multinomial log likelihood calculation for length frequency data
# with optional size selectivity

# Helper function to create test parameters
create_test_pars <- function() {
  list(
    k = 0.3,
    L_inf = 100,
    d = 0.2,
    m = 20,
    r = 0.5
  )
}

# Helper function to create test parameters with selectivity
create_test_pars_with_selectivity <- function() {
  list(
    k = 0.3,
    L_inf = 100,
    d = 0.2,
    m = 20,
    r = 0.5,
    l50 = 40,
    l25 = 30
  )
}

# Helper function to create test length frequency data
create_test_length_freq <- function() {
  data.frame(
    length = c(10, 20, 30, 40, 50, 60, 70),
    count = c(5, 15, 25, 30, 20, 10, 5)
  )
}

test_that("get_length_log_likelihood returns correct structure", {
  pars <- create_test_pars()
  length_freq <- create_test_length_freq()

  result <- get_length_log_likelihood(pars, length_freq)

  # Check that result is a single numeric value
  expect_type(result, "double")
  expect_length(result, 1)

  # Check that result is finite
  expect_true(is.finite(result))

  # Check that result is negative (log likelihood should be negative)
  expect_true(result < 0)
})

test_that("get_length_log_likelihood works without selectivity parameters", {
  pars <- create_test_pars()
  length_freq <- create_test_length_freq()

  # Should work fine without l50 and l25
  result <- get_length_log_likelihood(pars, length_freq)

  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("get_length_log_likelihood applies selectivity when parameters provided", {
  pars_no_sel <- create_test_pars()
  pars_with_sel <- create_test_pars_with_selectivity()
  length_freq <- create_test_length_freq()

  result_no_sel <- get_length_log_likelihood(pars_no_sel, length_freq)
  result_with_sel <- get_length_log_likelihood(pars_with_sel, length_freq)

  # Results should be different when selectivity is applied
  expect_false(result_no_sel == result_with_sel)

  # Both should be finite
  expect_true(is.finite(result_no_sel))
  expect_true(is.finite(result_with_sel))
})

test_that("get_length_log_likelihood selectivity function is correct", {
  # Test that selectivity values are correct at l25, l50, and l75
  pars <- create_test_pars_with_selectivity()

  # Calculate slope
  slope <- log(3) / (pars$l50 - pars$l25)

  # Test selectivity at key points
  sel_at_l25 <- 1 / (1 + exp(-slope * (pars$l25 - pars$l50)))
  sel_at_l50 <- 1 / (1 + exp(-slope * (pars$l50 - pars$l50)))
  sel_at_l75 <- 1 / (1 + exp(-slope * (pars$l50 + (pars$l50 - pars$l25) - pars$l50)))

  expect_equal(sel_at_l25, 0.25, tolerance = 1e-10)
  expect_equal(sel_at_l50, 0.50, tolerance = 1e-10)
  expect_equal(sel_at_l75, 0.75, tolerance = 1e-10)
})

test_that("get_length_log_likelihood skips selectivity with only l50", {
  pars_only_l50 <- create_test_pars()
  pars_only_l50$l50 <- 40

  pars_no_sel <- create_test_pars()
  length_freq <- create_test_length_freq()

  result_only_l50 <- get_length_log_likelihood(pars_only_l50, length_freq)
  result_no_sel <- get_length_log_likelihood(pars_no_sel, length_freq)

  # Should be identical when only l50 is provided (no selectivity applied)
  expect_equal(result_only_l50, result_no_sel)
})

test_that("get_length_log_likelihood skips selectivity with only l25", {
  pars_only_l25 <- create_test_pars()
  pars_only_l25$l25 <- 30

  pars_no_sel <- create_test_pars()
  length_freq <- create_test_length_freq()

  result_only_l25 <- get_length_log_likelihood(pars_only_l25, length_freq)
  result_no_sel <- get_length_log_likelihood(pars_no_sel, length_freq)

  # Should be identical when only l25 is provided (no selectivity applied)
  expect_equal(result_only_l25, result_no_sel)
})

test_that("get_length_log_likelihood handles different Delta_l values", {
  pars <- create_test_pars_with_selectivity()
  length_freq <- create_test_length_freq()

  result_1 <- get_length_log_likelihood(pars, length_freq, Delta_l = 1)
  result_2 <- get_length_log_likelihood(pars, length_freq, Delta_l = 2)
  result_5 <- get_length_log_likelihood(pars, length_freq, Delta_l = 5)

  # All should be finite
  expect_true(is.finite(result_1))
  expect_true(is.finite(result_2))
  expect_true(is.finite(result_5))

  # Results should be different (different discretization)
  expect_false(result_1 == result_2)
  expect_false(result_2 == result_5)
})

test_that("get_length_log_likelihood handles different l_max values", {
  pars <- create_test_pars_with_selectivity()
  length_freq <- create_test_length_freq()

  result_default <- get_length_log_likelihood(pars, length_freq)
  result_100 <- get_length_log_likelihood(pars, length_freq, l_max = 100)
  result_150 <- get_length_log_likelihood(pars, length_freq, l_max = 150)

  # All should be finite
  expect_true(is.finite(result_default))
  expect_true(is.finite(result_100))
  expect_true(is.finite(result_150))
})

test_that("get_length_log_likelihood handles different selectivity parameters", {
  length_freq <- create_test_length_freq()

  # Test different selectivity scenarios
  pars1 <- create_test_pars()
  pars1$l50 <- 30
  pars1$l25 <- 20

  pars2 <- create_test_pars()
  pars2$l50 <- 50
  pars2$l25 <- 40

  pars3 <- create_test_pars()
  pars3$l50 <- 40
  pars3$l25 <- 38  # Very steep

  result1 <- get_length_log_likelihood(pars1, length_freq)
  result2 <- get_length_log_likelihood(pars2, length_freq)
  result3 <- get_length_log_likelihood(pars3, length_freq)

  # All should be finite
  expect_true(is.finite(result1))
  expect_true(is.finite(result2))
  expect_true(is.finite(result3))

  # Results should be different
  expect_false(result1 == result2)
  expect_false(result2 == result3)
  expect_false(result1 == result3)
})

test_that("get_length_log_likelihood is deterministic", {
  pars <- create_test_pars_with_selectivity()
  length_freq <- create_test_length_freq()

  # Run multiple times with same parameters
  result1 <- get_length_log_likelihood(pars, length_freq)
  result2 <- get_length_log_likelihood(pars, length_freq)
  result3 <- get_length_log_likelihood(pars, length_freq)

  # Results should be identical
  expect_equal(result1, result2)
  expect_equal(result2, result3)
  expect_equal(result1, result3)
})

test_that("get_length_log_likelihood handles single observation", {
  pars <- create_test_pars_with_selectivity()
  length_freq <- data.frame(length = 30, count = 10)

  result <- get_length_log_likelihood(pars, length_freq)

  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("get_length_log_likelihood handles many observations", {
  pars <- create_test_pars_with_selectivity()

  # Create a larger dataset
  length_freq <- data.frame(
    length = 10:80,
    count = rpois(71, lambda = 20)
  )

  result <- get_length_log_likelihood(pars, length_freq)

  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("get_length_log_likelihood handles zero counts", {
  pars <- create_test_pars_with_selectivity()
  length_freq <- data.frame(
    length = c(20, 30, 40, 50),
    count = c(10, 0, 5, 0)
  )

  result <- get_length_log_likelihood(pars, length_freq)

  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("get_length_log_likelihood works with package data", {
  data(Cod_CS_pars)
  data(Cod_CS_length_freq)

  # Test without selectivity
  result_no_sel <- get_length_log_likelihood(Cod_CS_pars, Cod_CS_length_freq)

  expect_type(result_no_sel, "double")
  expect_true(is.finite(result_no_sel))
  expect_true(result_no_sel < 0)

  # Test with selectivity
  pars_with_sel <- Cod_CS_pars
  pars_with_sel$l50 <- 40
  pars_with_sel$l25 <- 30

  result_with_sel <- get_length_log_likelihood(pars_with_sel, Cod_CS_length_freq)

  expect_type(result_with_sel, "double")
  expect_true(is.finite(result_with_sel))
  expect_true(result_with_sel < 0)

  # Results should be different
  expect_false(result_no_sel == result_with_sel)
})

test_that("get_length_log_likelihood handles different parameter values", {
  length_freq <- create_test_length_freq()

  # Test with various parameter combinations
  test_pars_list <- list(
    list(k = 0.2, L_inf = 80, d = 0.1, m = 15, r = 0.3),
    list(k = 0.5, L_inf = 120, d = 0.3, m = 25, r = 0.6),
    list(k = 0.1, L_inf = 60, d = 0.05, m = 10, r = 0.2),
    list(k = 0.8, L_inf = 150, d = 0.4, m = 30, r = 0.8)
  )

  for (i in seq_along(test_pars_list)) {
    pars <- test_pars_list[[i]]
    result <- get_length_log_likelihood(pars, length_freq)

    expect_type(result, "double")
    expect_true(is.finite(result),
                info = paste("Non-finite result for parameter set", i))
  }
})

test_that("get_length_log_likelihood is sensitive to parameter changes", {
  pars_base <- create_test_pars_with_selectivity()
  length_freq <- create_test_length_freq()

  result_base <- get_length_log_likelihood(pars_base, length_freq)

  # Change k
  pars_k <- pars_base
  pars_k$k <- pars_k$k * 1.5
  result_k <- get_length_log_likelihood(pars_k, length_freq)
  expect_false(result_base == result_k)

  # Change L_inf
  pars_Linf <- pars_base
  pars_Linf$L_inf <- pars_Linf$L_inf * 1.2
  result_Linf <- get_length_log_likelihood(pars_Linf, length_freq)
  expect_false(result_base == result_Linf)

  # Change selectivity
  pars_sel <- pars_base
  pars_sel$l50 <- 50
  result_sel <- get_length_log_likelihood(pars_sel, length_freq)
  expect_false(result_base == result_sel)
})

test_that("get_length_log_likelihood validates selectivity increases with length", {
  # Manually verify selectivity function is increasing
  pars <- create_test_pars_with_selectivity()

  slope <- log(3) / (pars$l50 - pars$l25)
  lengths <- seq(10, 80, by = 10)

  selectivities <- numeric(length(lengths))
  for (i in seq_along(lengths)) {
    selectivities[i] <- 1 / (1 + exp(-slope * (lengths[i] - pars$l50)))
  }

  # Check that selectivity is increasing
  expect_true(all(diff(selectivities) > 0),
              "Selectivity should be strictly increasing with length")

  # Check bounds
  expect_true(all(selectivities >= 0 & selectivities <= 1),
              "Selectivity should be between 0 and 1")
})

test_that("get_length_log_likelihood handles extreme selectivity values", {
  length_freq <- create_test_length_freq()

  # Very steep selectivity
  pars_steep <- create_test_pars()
  pars_steep$l50 <- 40
  pars_steep$l25 <- 39.9
  result_steep <- get_length_log_likelihood(pars_steep, length_freq)
  expect_true(is.finite(result_steep))

  # Very gradual selectivity
  pars_gradual <- create_test_pars()
  pars_gradual$l50 <- 40
  pars_gradual$l25 <- 10
  result_gradual <- get_length_log_likelihood(pars_gradual, length_freq)
  expect_true(is.finite(result_gradual))
})

test_that("get_length_log_likelihood handles l_max smaller than max length", {
  pars <- create_test_pars_with_selectivity()
  length_freq <- data.frame(
    length = c(10, 20, 30, 40, 50),
    count = c(5, 10, 15, 10, 5)
  )

  # l_max smaller than max observed length (50)
  result <- get_length_log_likelihood(pars, length_freq, l_max = 40)

  # Should still work (max observed length gets mapped to grid)
  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("get_length_log_likelihood automatic l_max calculation", {
  pars <- create_test_pars()
  length_freq <- data.frame(
    length = c(10, 20, 30),
    count = c(5, 10, 5)
  )

  # Should automatically set l_max to 110% of max length (33)
  result <- get_length_log_likelihood(pars, length_freq)

  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("get_length_log_likelihood handles identical with and without vB_min_size", {
  pars_no_vB <- create_test_pars_with_selectivity()
  pars_with_vB <- create_test_pars_with_selectivity()
  pars_with_vB$vB_min_size <- 10

  length_freq <- create_test_length_freq()

  result_no_vB <- get_length_log_likelihood(pars_no_vB, length_freq)
  result_with_vB <- get_length_log_likelihood(pars_with_vB, length_freq)

  # Both should be finite
  expect_true(is.finite(result_no_vB))
  expect_true(is.finite(result_with_vB))

  # Results should be different (different growth model)
  expect_false(result_no_vB == result_with_vB)
})

