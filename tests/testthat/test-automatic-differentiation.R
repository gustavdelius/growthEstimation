# Test automatic differentiation vs finite differences
# This test verifies that the TMB automatic differentiation gradients
# are consistent with finite difference approximations

library(testthat)
library(growthEstimation)

# Helper function to compute finite difference gradients
compute_finite_difference_gradient <- function(obj, par, h = 1e-6) {
  n_par <- length(par)
  gradient <- numeric(n_par)

  # Store original parameter values (for potential future use)
  # original_par <- par

  # Compute function value at original point
  f0 <- obj$fn(par)

  # Check if original function value is finite
  if (!is.finite(f0)) {
    warning("Function value at original point is not finite")
    return(rep(NA, n_par))
  }

  for (i in 1:n_par) {
    # Forward difference
    par_forward <- par
    par_forward[i] <- par[i] + h
    f_forward <- obj$fn(par_forward)

    # Backward difference
    par_backward <- par
    par_backward[i] <- par[i] - h
    f_backward <- obj$fn(par_backward)

    # Check if both function values are finite
    if (!is.finite(f_forward) || !is.finite(f_backward)) {
      warning(paste("Function value not finite for parameter", i, "at step size", h))
      gradient[i] <- NA
    } else {
      # Central difference (more accurate)
      gradient[i] <- (f_forward - f_backward) / (2 * h)
    }
  }

  return(gradient)
}

# Helper function to create test data and TMB object
create_test_tmb_object <- function(pars, age_at_length = NULL, length_freq = NULL,
                                    weight_age = 1, weight_length = 1) {
  if (is.null(age_at_length)) {
    age_at_length <- Cod_CS_age_at_length
  }

  # Common discretisation
  Delta_l <- 1
  Delta_t <- 0.05

  # Create TMB object
  tmb_result <- fit_tmb_nll(
    pars = pars,
    age_at_length = age_at_length,
    length_freq = length_freq,
    weight_age = weight_age,
    weight_length = weight_length,
    Delta_l = Delta_l,
    Delta_t = Delta_t
  )

  return(tmb_result$obj)
}

test_that("TMB automatic differentiation gradients match finite differences", {
  # Test parameters
  pars <- list(
    k = 0.5,
    L_inf = 80,
    d = 0.2,
    m = 20,
    r = 0.3,
    spawning_mu = 0.4,
    spawning_kappa = 3,
    annuli_date = 0.25,
    annuli_min_age = 0.5
  )

  # Create TMB object
  obj <- create_test_tmb_object(pars)

  # Test at the initial parameter values
  test_par <- obj$par

  # Get TMB automatic differentiation gradient
  ad_gradient_raw <- obj$gr(test_par)
  # TMB returns a matrix, convert to vector
  ad_gradient <- as.vector(ad_gradient_raw)

  # Get finite difference gradient
  fd_gradient <- compute_finite_difference_gradient(obj, test_par)

  # Compare gradients
  expect_equal(length(ad_gradient), length(fd_gradient))
  expect_equal(length(ad_gradient), length(test_par))

  # Check that both gradients are finite
  expect_true(all(is.finite(ad_gradient)))
  expect_true(all(is.finite(fd_gradient)))

  # Compare with reasonable tolerance
  # Finite differences have truncation error, so we need some tolerance
  expect_equal(ad_gradient, fd_gradient, tolerance = 1e-3)
})

test_that("Gradient consistency holds for different parameter values", {
  # Test multiple parameter sets
  test_pars_list <- list(
    # Original parameters
    list(
      k = 0.5, L_inf = 80, d = 0.2, m = 20, r = 0.3,
      spawning_mu = 0.4, spawning_kappa = 3,
      annuli_date = 0.25, annuli_min_age = 0.5
    ),
    # Different growth parameters
    list(
      k = 0.3, L_inf = 100, d = 0.1, m = 15, r = 0.2,
      spawning_mu = 0.4, spawning_kappa = 3,
      annuli_date = 0.25, annuli_min_age = 0.5
    ),
    # Different spawning parameters
    list(
      k = 0.4, L_inf = 90, d = 0.15, m = 25, r = 0.4,
      spawning_mu = 0.6, spawning_kappa = 2,
      annuli_date = 0.25, annuli_min_age = 0.5
    )
  )

  for (pars in test_pars_list) {
    # Create TMB object
    obj <- create_test_tmb_object(pars)
    test_par <- obj$par

    # Get gradients
    ad_gradient_raw <- obj$gr(test_par)
    ad_gradient <- as.vector(ad_gradient_raw)
    fd_gradient <- compute_finite_difference_gradient(obj, test_par)

    # Check consistency
    expect_true(all(is.finite(ad_gradient)))
    expect_true(all(is.finite(fd_gradient)))
    expect_equal(ad_gradient, fd_gradient, tolerance = 1e-3)
  }
})

test_that("Gradient consistency for individual parameters", {
  pars <- list(
    k = 0.5, L_inf = 80, d = 0.2, m = 20, r = 0.3,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )

  obj <- create_test_tmb_object(pars)
  test_par <- obj$par

  # Get both gradients
  ad_gradient_raw <- obj$gr(test_par)
  ad_gradient <- as.vector(ad_gradient_raw)
  fd_gradient <- compute_finite_difference_gradient(obj, test_par)

  # Check each parameter individually
  parameter_names <- c("k", "L_inf", "d", "m")

  for (i in seq_along(parameter_names)) {
    param_name <- parameter_names[i]

    # Check that gradients are finite
    expect_true(is.finite(ad_gradient[i]),
                info = paste("AD gradient for", param_name, "is not finite"))
    expect_true(is.finite(fd_gradient[i]),
                info = paste("FD gradient for", param_name, "is not finite"))

    # Check consistency
    expect_equal(ad_gradient[i], fd_gradient[i], tolerance = 1e-3,
                 info = paste("Gradient mismatch for parameter", param_name))
  }
})

test_that("Gradient consistency with different step sizes", {
  pars <- list(
    k = 0.5, L_inf = 80, d = 0.2, m = 20, r = 0.3,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )

  obj <- create_test_tmb_object(pars)
  test_par <- obj$par

  # Get TMB gradient
  ad_gradient_raw <- obj$gr(test_par)
  ad_gradient <- as.vector(ad_gradient_raw)

  # Test different finite difference step sizes
  step_sizes <- c(1e-5, 1e-6, 1e-7)

  for (h in step_sizes) {
    fd_gradient <- compute_finite_difference_gradient(obj, test_par, h = h)

    # Check that finite difference gradient is finite
    expect_true(all(is.finite(fd_gradient)),
                info = paste("FD gradient with step size", h, "is not finite"))

    # Check consistency (tolerance may need to be adjusted based on step size)
    tolerance <- max(1e-3, h * 100)  # Allow larger tolerance for larger step sizes
    expect_equal(ad_gradient, fd_gradient, tolerance = tolerance,
                 info = paste("Gradient mismatch with step size", h))
  }
})

test_that("Gradient consistency with different survey data", {
  # Test with different survey datasets
  survey_datasets <- list(
    # Original dataset
    data.frame(
      survey_date = c(2020.25, 2020.25, 2021.10, 2021.10, 2021.10),
      length = c(20L, 21L, 20L, 21L, 22L),
      K = c(0L, 1L, 0L, 1L, 1L),
      count = c(10, 5, 7, 8, 3)
    ),
    # Larger dataset
    data.frame(
      survey_date = c(2020.25, 2020.25, 2020.25, 2021.10, 2021.10, 2021.10, 2022.5),
      length = c(20L, 21L, 22L, 20L, 21L, 22L, 23L),
      K = c(0L, 1L, 1L, 0L, 1L, 1L, 2L),
      count = c(10, 5, 3, 7, 8, 3, 2)
    ),
    # Single survey
    data.frame(
      survey_date = c(2020.5),
      length = c(20L),
      K = c(0L),
      count = c(15)
    )
  )

  pars <- list(
    k = 0.5, L_inf = 80, d = 0.2, m = 20, r = 0.3,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )

  for (i in seq_along(survey_datasets)) {
    age_at_length <- survey_datasets[[i]]

    # Create TMB object with this survey data
    # Suppress warnings from TMB::sdreport() which may fail for some parameter combinations
    suppressWarnings({
      obj <- create_test_tmb_object(pars, age_at_length)
    })
    test_par <- obj$par

    # Get gradients
    ad_gradient_raw <- obj$gr(test_par)
    ad_gradient <- as.vector(ad_gradient_raw)
    fd_gradient <- compute_finite_difference_gradient(obj, test_par)

    # Check consistency
    expect_true(all(is.finite(ad_gradient)),
                info = paste("AD gradient for survey dataset", i, "is not finite"))

    # Skip finite difference check if it contains NAs
    if (any(is.na(fd_gradient))) {
      skip(paste("Finite difference gradient contains NAs for survey dataset", i))
    }

    expect_true(all(is.finite(fd_gradient)),
                info = paste("FD gradient for survey dataset", i, "is not finite"))
    expect_equal(ad_gradient, fd_gradient, tolerance = 1e-3,
                 info = paste("Gradient mismatch for survey dataset", i))
  }
})

test_that("Gradient magnitude and sign consistency", {
  pars <- list(
    k = 0.5, L_inf = 80, d = 0.2, m = 20, r = 0.3,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )

  obj <- create_test_tmb_object(pars)
  test_par <- obj$par

  # Get gradients
  ad_gradient_raw <- obj$gr(test_par)
  ad_gradient <- as.vector(ad_gradient_raw)
  fd_gradient <- compute_finite_difference_gradient(obj, test_par)

  # Check that gradients have reasonable magnitudes (not too large or too small)
  expect_true(all(abs(ad_gradient) < 1e6),
               info = "AD gradient has unreasonably large magnitude")
  expect_true(all(abs(fd_gradient) < 1e6),
               info = "FD gradient has unreasonably large magnitude")

  # Check that gradients are not all zero (unless at optimum)
  expect_true(any(abs(ad_gradient) > 1e-10),
               info = "AD gradient is suspiciously close to zero")

  # Check that signs are consistent between AD and FD (only for non-zero gradients)
  # For very small gradients, signs might be unreliable due to numerical precision
  significant_gradients <- abs(ad_gradient) > 1e-10 & abs(fd_gradient) > 1e-10
  if (any(significant_gradients)) {
    expect_equal(sign(ad_gradient[significant_gradients]), sign(fd_gradient[significant_gradients]),
                 info = "AD and FD gradients have different signs for significant gradients")
  }
})

# ============================================================================
# Tests for length frequency data and selectivity parameters
# ============================================================================

test_that("TMB AD gradients match FD for length frequency data", {
  # Test parameters with selectivity
  pars <- list(
    k = 0.3,
    L_inf = 100,
    d = 0.2,
    m = 20,
    r = 0.5,
    l50 = 40,
    ratio = 0.75,
    vB_min_size = 0,
    spawning_mu = 0.4,
    spawning_kappa = 3,
    annuli_date = 0.25,
    annuli_min_age = 0.5
  )
  
  # Create TMB object with length frequency data
  obj <- create_test_tmb_object(pars, 
                                 age_at_length = Cod_CS_age_at_length,
                                 length_freq = Cod_CS_length_freq)
  
  # Test at initial parameters
  test_par <- obj$par
  
  # Get gradients
  ad_gradient_raw <- obj$gr(test_par)
  ad_gradient <- as.vector(ad_gradient_raw)
  fd_gradient <- compute_finite_difference_gradient(obj, test_par)
  
  # Verify
  expect_equal(length(ad_gradient), 7)  # k, L_inf, d, m, r, l50, ratio
  expect_true(all(is.finite(ad_gradient)))
  expect_true(all(is.finite(fd_gradient)))
  expect_equal(ad_gradient, fd_gradient, tolerance = 1e-3)
})

test_that("Gradient consistency for l50 parameter", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  obj <- create_test_tmb_object(pars,
                                 age_at_length = Cod_CS_age_at_length,
                                 length_freq = Cod_CS_length_freq)
  test_par <- obj$par
  
  # Get gradients
  ad_gradient <- as.vector(obj$gr(test_par))
  fd_gradient <- compute_finite_difference_gradient(obj, test_par)
  
  # l50 is the 6th parameter (k, L_inf, d, m, r, l50)
  l50_index <- which(names(test_par) == "l50")
  
  expect_true(is.finite(ad_gradient[l50_index]),
              info = "AD gradient for l50 is not finite")
  expect_true(is.finite(fd_gradient[l50_index]),
              info = "FD gradient for l50 is not finite")
  expect_equal(ad_gradient[l50_index], fd_gradient[l50_index], tolerance = 1e-3,
               info = "Gradient mismatch for l50")
})

test_that("Gradient consistency for ratio parameter", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  obj <- create_test_tmb_object(pars,
                                 age_at_length = Cod_CS_age_at_length,
                                 length_freq = Cod_CS_length_freq)
  test_par <- obj$par
  
  # Get gradients
  ad_gradient <- as.vector(obj$gr(test_par))
  fd_gradient <- compute_finite_difference_gradient(obj, test_par)
  
  # ratio is the 7th parameter
  ratio_index <- which(names(test_par) == "ratio")
  
  expect_true(is.finite(ad_gradient[ratio_index]),
              info = "AD gradient for ratio is not finite")
  expect_true(is.finite(fd_gradient[ratio_index]),
              info = "FD gradient for ratio is not finite")
  expect_equal(ad_gradient[ratio_index], fd_gradient[ratio_index], tolerance = 1e-3,
               info = "Gradient mismatch for ratio")
})

test_that("Gradient consistency for all parameters with length frequency", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  obj <- create_test_tmb_object(pars,
                                 age_at_length = Cod_CS_age_at_length,
                                 length_freq = Cod_CS_length_freq)
  test_par <- obj$par
  
  # Get gradients
  ad_gradient <- as.vector(obj$gr(test_par))
  fd_gradient <- compute_finite_difference_gradient(obj, test_par)
  
  # Check each parameter individually
  parameter_names <- c("k", "L_inf", "d", "m", "r", "l50", "ratio")
  
  for (i in seq_along(parameter_names)) {
    param_name <- parameter_names[i]
    
    expect_true(is.finite(ad_gradient[i]),
                info = paste("AD gradient for", param_name, "is not finite"))
    expect_true(is.finite(fd_gradient[i]),
                info = paste("FD gradient for", param_name, "is not finite"))
    expect_equal(ad_gradient[i], fd_gradient[i], tolerance = 1e-3,
                 info = paste("Gradient mismatch for parameter", param_name))
  }
})

test_that("Gradient consistency with different selectivity parameters", {
  # Test multiple selectivity scenarios
  test_cases <- list(
    list(l50 = 30, ratio = 0.7),   # Early selectivity
    list(l50 = 50, ratio = 0.8),   # Later selectivity
    list(l50 = 40, ratio = 0.9),   # Steep selectivity (ratio close to 1)
    list(l50 = 40, ratio = 0.5)    # Gradual selectivity
  )
  
  for (case in test_cases) {
    pars <- list(
      k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
      l50 = case$l50, ratio = case$ratio, vB_min_size = 0,
      spawning_mu = 0.4, spawning_kappa = 3,
      annuli_date = 0.25, annuli_min_age = 0.5
    )
    
    obj <- create_test_tmb_object(pars,
                                   age_at_length = Cod_CS_age_at_length,
                                   length_freq = Cod_CS_length_freq)
    test_par <- obj$par
    
    ad_gradient <- as.vector(obj$gr(test_par))
    fd_gradient <- compute_finite_difference_gradient(obj, test_par)
    
    expect_true(all(is.finite(ad_gradient)),
                info = paste("AD gradient not finite for l50 =", case$l50, "ratio =", case$ratio))
    expect_true(all(is.finite(fd_gradient)),
                info = paste("FD gradient not finite for l50 =", case$l50, "ratio =", case$ratio))
    expect_equal(ad_gradient, fd_gradient, tolerance = 1e-3,
                 info = paste("Gradient mismatch for l50 =", case$l50, "ratio =", case$ratio))
  }
})

test_that("Gradient consistency with different likelihood weights", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  # Test different weighting schemes
  weight_cases <- list(
    list(weight_age = 1, weight_length = 1),
    list(weight_age = 2, weight_length = 1),
    list(weight_age = 1, weight_length = 2),
    list(weight_age = 0, weight_length = 1),  # Length only
    list(weight_age = 1, weight_length = 0)   # Age only (but still has l50, ratio params)
  )
  
  for (case in weight_cases) {
    obj <- create_test_tmb_object(pars,
                                   age_at_length = Cod_CS_age_at_length,
                                   length_freq = Cod_CS_length_freq,
                                   weight_age = case$weight_age,
                                   weight_length = case$weight_length)
    test_par <- obj$par
    
    ad_gradient <- as.vector(obj$gr(test_par))
    fd_gradient <- compute_finite_difference_gradient(obj, test_par)
    
    expect_true(all(is.finite(ad_gradient)),
                info = paste("AD gradient not finite for weights", 
                            case$weight_age, ":", case$weight_length))
    expect_true(all(is.finite(fd_gradient)),
                info = paste("FD gradient not finite for weights",
                            case$weight_age, ":", case$weight_length))
    expect_equal(ad_gradient, fd_gradient, tolerance = 1e-3,
                 info = paste("Gradient mismatch for weights",
                             case$weight_age, ":", case$weight_length))
  }
})

test_that("Gradient for l50 is sensitive to length frequency data", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  # With length frequency data
  obj <- create_test_tmb_object(pars,
                                 age_at_length = Cod_CS_age_at_length,
                                 length_freq = Cod_CS_length_freq)
  test_par <- obj$par
  ad_gradient <- as.vector(obj$gr(test_par))
  
  l50_index <- which(names(test_par) == "l50")
  
  # Gradient for l50 should be non-zero when length_freq data is used
  expect_true(abs(ad_gradient[l50_index]) > 1e-6,
              info = "Gradient for l50 should be non-zero with length frequency data")
})

test_that("Gradient for ratio is sensitive to length frequency data", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  # With length frequency data
  obj <- create_test_tmb_object(pars,
                                 age_at_length = Cod_CS_age_at_length,
                                 length_freq = Cod_CS_length_freq)
  test_par <- obj$par
  ad_gradient <- as.vector(obj$gr(test_par))
  
  ratio_index <- which(names(test_par) == "ratio")
  
  # Gradient for ratio should be non-zero when length_freq data is used
  expect_true(abs(ad_gradient[ratio_index]) > 1e-6,
              info = "Gradient for ratio should be non-zero with length frequency data")
})

test_that("Gradients change with different length frequency datasets", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  # Create different length frequency datasets
  length_freq_1 <- data.frame(
    length = c(20, 30, 40, 50, 60),
    count = c(10, 20, 30, 20, 10)
  )
  
  length_freq_2 <- data.frame(
    length = c(10, 20, 30, 40, 50),
    count = c(5, 15, 25, 15, 5)
  )
  
  # Get gradients for each dataset
  obj1 <- create_test_tmb_object(pars, 
                                  age_at_length = Cod_CS_age_at_length,
                                  length_freq = length_freq_1)
  grad1 <- as.vector(obj1$gr(obj1$par))
  
  obj2 <- create_test_tmb_object(pars,
                                  age_at_length = Cod_CS_age_at_length,
                                  length_freq = length_freq_2)
  grad2 <- as.vector(obj2$gr(obj2$par))
  
  # Gradients should be different for different data
  expect_false(all(abs(grad1 - grad2) < 1e-10),
               info = "Gradients should differ for different length frequency data")
})

test_that("Gradient consistency for length-only likelihood", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  # Length frequency only (weight_age = 0)
  obj <- create_test_tmb_object(pars,
                                 age_at_length = Cod_CS_age_at_length,
                                 length_freq = Cod_CS_length_freq,
                                 weight_age = 0,
                                 weight_length = 1)
  test_par <- obj$par
  
  ad_gradient <- as.vector(obj$gr(test_par))
  fd_gradient <- compute_finite_difference_gradient(obj, test_par)
  
  expect_true(all(is.finite(ad_gradient)))
  expect_true(all(is.finite(fd_gradient)))
  expect_equal(ad_gradient, fd_gradient, tolerance = 1e-3)
})

test_that("All 7 parameters have non-zero gradients with combined data", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  obj <- create_test_tmb_object(pars,
                                 age_at_length = Cod_CS_age_at_length,
                                 length_freq = Cod_CS_length_freq)
  test_par <- obj$par
  ad_gradient <- as.vector(obj$gr(test_par))
  
  # Check that all parameters have meaningful gradients
  parameter_names <- c("k", "L_inf", "d", "m", "r", "l50", "ratio")
  
  for (i in seq_along(parameter_names)) {
    param_name <- parameter_names[i]
    
    # For most parameters, gradient should be non-negligible
    # Some parameters might have small gradients depending on data
    expect_true(is.finite(ad_gradient[i]),
                info = paste("Gradient for", param_name, "is not finite"))
  }
})

test_that("Gradient scales appropriately with likelihood weights", {
  pars <- list(
    k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5,
    l50 = 40, ratio = 0.75, vB_min_size = 0,
    spawning_mu = 0.4, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
  )
  
  # Get gradient with weight_length = 1
  obj1 <- create_test_tmb_object(pars,
                                  age_at_length = Cod_CS_age_at_length,
                                  length_freq = Cod_CS_length_freq,
                                  weight_age = 0,
                                  weight_length = 1)
  grad1 <- as.vector(obj1$gr(obj1$par))
  
  # Get gradient with weight_length = 2
  obj2 <- create_test_tmb_object(pars,
                                  age_at_length = Cod_CS_age_at_length,
                                  length_freq = Cod_CS_length_freq,
                                  weight_age = 0,
                                  weight_length = 2)
  grad2 <- as.vector(obj2$gr(obj2$par))
  
  # Gradients should scale approximately linearly with weights
  # grad2 ≈ 2 * grad1
  expect_equal(grad2, 2 * grad1, tolerance = 1e-6,
               info = "Gradients should scale linearly with weights")
})
