# Test suite for solve_pde_steady_state function
# Tests the steady state solution of the fish abundance PDE

library(testthat)
library(growthEstimation)

# Helper function to create test parameters
create_test_pars <- function() {
  list(
    k = 0.5,
    L_inf = 80,
    d = 0.2,
    m = 20,
    r = 0.3,
    vB_min_size = 0
  )
}

# Helper function to create test parameters with different vB_min_size
create_test_pars_vB <- function() {
  list(
    k = 0.4,
    L_inf = 100,
    d = 0.15,
    m = 25,
    r = 0.2,
    vB_min_size = 10
  )
}

test_that("solve_pde_steady_state returns correct structure", {
  pars <- create_test_pars()
  result <- solve_pde_steady_state(pars, Delta_l = 2, l_max = 50)
  
  # Check that result is a numeric vector
  expect_type(result, "double")
  expect_true(is.vector(result))
  
  # Check that result has correct length (should be N_l = ceiling(l_max / Delta_l))
  expected_length <- ceiling(50 / 2)  # 25
  expect_equal(length(result), expected_length)
  
  # Check that all values are finite
  expect_true(all(is.finite(result)))
  
  # Check that result is normalized (max absolute value should be 1)
  expect_equal(max(abs(result)), 1.0)
})

test_that("solve_pde_steady_state handles different grid sizes", {
  pars <- create_test_pars()
  
  # Test different combinations of Delta_l and l_max
  test_cases <- list(
    list(Delta_l = 1, l_max = 20),
    list(Delta_l = 2, l_max = 40),
    list(Delta_l = 5, l_max = 50),
    list(Delta_l = 0.5, l_max = 10)
  )
  
  for (case in test_cases) {
    result <- solve_pde_steady_state(pars, Delta_l = case$Delta_l, l_max = case$l_max)
    
    expected_length <- ceiling(case$l_max / case$Delta_l)
    expect_equal(length(result), expected_length,
                 info = paste("Failed for Delta_l =", case$Delta_l, "l_max =", case$l_max))
    
    expect_true(all(is.finite(result)),
                info = paste("Non-finite values for Delta_l =", case$Delta_l, "l_max =", case$l_max))
    
    expect_equal(max(abs(result)), 1.0,
                 info = paste("Not normalized for Delta_l =", case$Delta_l, "l_max =", case$l_max))
  }
})

test_that("solve_pde_steady_state handles different parameter sets", {
  # Test with different parameter combinations
  test_pars_list <- list(
    create_test_pars(),
    create_test_pars_vB(),
    list(k = 0.3, L_inf = 60, d = 0.1, m = 15, r = 0.4, vB_min_size = 5),
    list(k = 0.8, L_inf = 120, d = 0.3, m = 30, r = 0.1, vB_min_size = 15)
  )
  
  for (i in seq_along(test_pars_list)) {
    pars <- test_pars_list[[i]]
    result <- solve_pde_steady_state(pars, Delta_l = 2, l_max = 50)
    
    expect_type(result, "double")
    expect_true(all(is.finite(result)),
                info = paste("Non-finite values for parameter set", i))
    expect_equal(max(abs(result)), 1.0,
                 info = paste("Not normalized for parameter set", i))
  }
})

test_that("solve_pde_steady_state handles boundary conditions correctly", {
  pars <- create_test_pars()
  result <- solve_pde_steady_state(pars, Delta_l = 1, l_max = 20)
  
  # Check that the solution doesn't have extreme values at boundaries
  # (exact values depend on the specific PDE, but should be reasonable)
  expect_true(abs(result[1]) <= 1.0, "First value should be within normalized range")
  expect_true(abs(result[length(result)]) <= 1.0, "Last value should be within normalized range")
  
  # Check that solution is not identically zero
  expect_true(max(abs(result)) > 1e-10, "Solution should not be identically zero")
})

test_that("solve_pde_steady_state handles edge cases", {
  pars <- create_test_pars()
  
  # Test with very small domain
  result_small <- solve_pde_steady_state(pars, Delta_l = 1, l_max = 1)
  expect_equal(length(result_small), 1)
  expect_true(is.finite(result_small))
  
  # Test with very large step size
  result_large_step <- solve_pde_steady_state(pars, Delta_l = 10, l_max = 50)
  expect_equal(length(result_large_step), 5)  # ceiling(50/10) = 5
  expect_true(all(is.finite(result_large_step)))
  
  # Test with parameters that might cause numerical issues
  pars_extreme <- list(k = 0.01, L_inf = 200, d = 0.001, m = 100, r = 0.01, vB_min_size = 0)
  result_extreme <- solve_pde_steady_state(pars_extreme, Delta_l = 2, l_max = 100)
  expect_true(all(is.finite(result_extreme)))
  expect_equal(max(abs(result_extreme)), 1.0)
})

test_that("solve_pde_steady_state is consistent with different vB_min_size values", {
  pars_base <- create_test_pars()
  
  # Test with vB_min_size = 0 (all von Bertalanffy)
  pars_vB0 <- pars_base
  pars_vB0$vB_min_size <- 0
  result_vB0 <- solve_pde_steady_state(pars_vB0, Delta_l = 2, l_max = 50)
  
  # Test with vB_min_size = 20 (mixed growth)
  pars_vB20 <- pars_base
  pars_vB20$vB_min_size <- 20
  result_vB20 <- solve_pde_steady_state(pars_vB20, Delta_l = 2, l_max = 50)
  
  # Test with vB_min_size = 50 (all constant growth)
  pars_vB50 <- pars_base
  pars_vB50$vB_min_size <- 50
  result_vB50 <- solve_pde_steady_state(pars_vB50, Delta_l = 2, l_max = 50)
  
  # All should produce valid solutions
  expect_true(all(is.finite(result_vB0)))
  expect_true(all(is.finite(result_vB20)))
  expect_true(all(is.finite(result_vB50)))
  
  # All should be normalized
  expect_equal(max(abs(result_vB0)), 1.0)
  expect_equal(max(abs(result_vB20)), 1.0)
  expect_equal(max(abs(result_vB50)), 1.0)
  
  # Solutions should be different (unless coincidentally the same)
  # We don't require them to be different, just that they're valid
})

test_that("solve_pde_steady_state handles missing vB_min_size parameter", {
  # Test with pars that don't have vB_min_size (should default to 0)
  pars_no_vB <- list(k = 0.5, L_inf = 80, d = 0.2, m = 20, r = 0.3)
  
  result <- solve_pde_steady_state(pars_no_vB, Delta_l = 2, l_max = 50)
  
  expect_type(result, "double")
  expect_true(all(is.finite(result)))
  expect_equal(max(abs(result)), 1.0)
})

test_that("solve_pde_steady_state produces reasonable solutions", {
  pars <- create_test_pars()
  result <- solve_pde_steady_state(pars, Delta_l = 1, l_max = 50)
  
  # Check that solution doesn't have unreasonable oscillations
  # (This is a qualitative check - the solution should be "smooth" in some sense)
  
  # Check that consecutive values don't differ by more than a reasonable amount
  diff_result <- abs(diff(result))
  expect_true(max(diff_result) < 10.0, "Solution should not have extreme jumps")
  
  # Check that solution doesn't have too many sign changes
  sign_changes <- sum(diff(sign(result)) != 0)
  expect_true(sign_changes < length(result) / 2, 
              "Solution should not oscillate excessively")
})

test_that("solve_pde_steady_state is deterministic", {
  pars <- create_test_pars()
  
  # Run the function multiple times with the same parameters
  result1 <- solve_pde_steady_state(pars, Delta_l = 2, l_max = 50)
  result2 <- solve_pde_steady_state(pars, Delta_l = 2, l_max = 50)
  result3 <- solve_pde_steady_state(pars, Delta_l = 2, l_max = 50)
  
  # Results should be identical (within numerical precision)
  expect_equal(result1, result2, tolerance = 1e-15)
  expect_equal(result2, result3, tolerance = 1e-15)
  expect_equal(result1, result3, tolerance = 1e-15)
})

test_that("solve_pde_steady_state handles zero parameters gracefully", {
  # Test with parameters that might cause division by zero
  pars_zero <- list(k = 0, L_inf = 80, d = 0, m = 0, r = 0.3, vB_min_size = 0)
  
  # This should not crash, but may produce a degenerate solution
  result <- solve_pde_steady_state(pars_zero, Delta_l = 2, l_max = 20)
  
  # Result should still be finite and normalized
  expect_true(all(is.finite(result)))
  expect_equal(max(abs(result)), 1.0)
})

test_that("solve_pde_steady_state works with package data", {
  # Test with the actual package data
  data(Cod_CS_pars)
  
  result <- solve_pde_steady_state(Cod_CS_pars, Delta_l = 2, l_max = 100)
  
  expect_type(result, "double")
  expect_true(all(is.finite(result)))
  expect_equal(max(abs(result)), 1.0)
  expect_equal(length(result), ceiling(100 / 2))  # Should be 50
})
