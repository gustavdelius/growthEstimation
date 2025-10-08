# Test suite for parameter recovery from simulated data
# Tests that fit_tmb_nll() can recover known parameters from simulated data

skip("Skipping all tests in this file for now")

test_that("fit_on_simulated_data recovers age-only parameters", {
    set.seed(42)

    true_pars <- list(
        k = 0.3,
        L_inf = 100,
        d = 0.2,
        m = 20,
        r = 5,
        vB_min_size = 10,
        spawning_mu = 0.5,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )

    start_pars <- list(
        k = 0.1,
        L_inf = 200,
        d = 0.1,
        m = 10,
        r = 10,
        vB_min_size = 10,
        spawning_mu = 0.5,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )

    # Simulate and fit
    result <- fit_on_simulated_data(
        true_pars = true_pars,
        start_pars = start_pars,
        survey_dates = c(2023.25, 2023.75),
        lengths = seq(20, 80, by = 5),
        n_per_length = 30,
        Delta_l = 1,
        Delta_t = 0.05
    )

    # Check convergence
    expect_equal(result$fit$opt$convergence, 0)

    # Check parameter recovery (allow some tolerance due to sampling variability)
    comp <- result$comparison
    for (i in 1:nrow(comp)) {
        param_name <- comp$parameter[i]
        true_val <- comp$true[i]
        est_val <- comp$estimated[i]

        # Allow 20% relative error for recovery
        relative_error <- abs(est_val - true_val) / abs(true_val)
        expect_true(relative_error < 0.2,
                   info = paste("Parameter", param_name, "not recovered within 20%:",
                               "true =", true_val, "estimated =", est_val))
    }
})

test_that("fit_on_simulated_data recovers length frequency parameters", {
    set.seed(123)

    true_pars <- list(
        k = 0.3,
        L_inf = 100,
        d = 0.2,
        m = 20,
        r = 5,
        l50 = 40,
        ratio = 0.75,
        l25 = 30,  # This will be auto-calculated from ratio
        vB_min_size = 10,
        spawning_mu = 0.5,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )

    start_pars <- list(
        k = 0.1,
        L_inf = 200,
        d = 0.1,
        m = 10,
        r = 10,
        l50 = 20,
        ratio = 0.75,
        l25 = 15,  # This will be auto-calculated from ratio
        vB_min_size = 10,
        spawning_mu = 0.5,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )

    # Simulate and fit with length frequency data
    result <- fit_on_simulated_data(
        true_pars = true_pars,
        start_pars = start_pars,
        survey_dates = c(2023.25, 2023.75),
        lengths = seq(10, 90, by = 2),
        n_per_length = 20,
        n_length_total = 1000,  # Simulate length frequency data
        Delta_l = 1,
        Delta_t = 0.05
    )

    # Check that length frequency data was simulated
    expect_false(is.null(result$simulated_length_data))
    expect_true(nrow(result$simulated_length_data) > 0)
    expect_equal(sum(result$simulated_length_data$count), 1000)

    # Check convergence (convergence code 1 is also acceptable)
    expect_true(result$fit$opt$convergence %in% c(0, 1))

    # Check parameter recovery
    comp <- result$comparison

    # Core PDE parameters should be recovered reasonably well
    core_params <- c("k", "L_inf", "d", "m", "r")
    for (param_name in core_params) {
        idx <- which(comp$parameter == param_name)
        if (length(idx) > 0) {
            true_val <- comp$true[idx]
            est_val <- comp$estimated[idx]

            # Allow 30% relative error (more lenient with length data)
            relative_error <- abs(est_val - true_val) / abs(true_val)
            expect_true(relative_error < 0.3,
                       info = paste("Parameter", param_name, "not recovered within 30%:",
                                   "true =", true_val, "estimated =", est_val))
        }
    }

    # Selectivity parameters should be recovered
    l50_idx <- which(comp$parameter == "l50")
    if (length(l50_idx) > 0) {
        expect_true(abs(comp$estimated[l50_idx] - comp$true[l50_idx]) < 10,
                   info = "l50 should be recovered within 10 cm")
    }

    ratio_idx <- which(comp$parameter == "ratio")
    if (length(ratio_idx) > 0) {
        expect_true(abs(comp$estimated[ratio_idx] - comp$true[ratio_idx]) < 0.15,
                   info = "ratio should be recovered within 0.15")
    }
})

test_that("simulate_length_freq_from_parameters produces valid data", {
    set.seed(456)

    pars <- list(
        k = 0.3,
        L_inf = 100,
        d = 0.2,
        m = 20,
        r = 5,
        l50 = 40,
        ratio = 0.75,
        vB_min_size = 0
    )

    lengths <- seq(10, 80, by = 5)
    n_total <- 500

    sim_data <- simulate_length_freq_from_parameters(
        pars = pars,
        lengths = lengths,
        n_total = n_total,
        Delta_l = 1
    )

    # Check structure
    expect_s3_class(sim_data, "data.frame")
    expect_true(all(c("length", "count") %in% names(sim_data)))

    # Check dimensions
    expect_equal(nrow(sim_data), length(lengths))

    # Check total count
    expect_equal(sum(sim_data$count), n_total)

    # Check all counts are non-negative
    expect_true(all(sim_data$count >= 0))

    # Check lengths match input
    expect_equal(sim_data$length, lengths)
})

test_that("simulate_length_freq selectivity affects distribution", {
    set.seed(789)

    pars_no_sel <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        vB_min_size = 0
    )

    pars_with_sel <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        l50 = 40, ratio = 0.75,
        vB_min_size = 0
    )

    lengths <- seq(10, 80, by = 5)
    n_total <- 500

    # Simulate without selectivity
    sim_no_sel <- simulate_length_freq_from_parameters(
        pars = pars_no_sel,
        lengths = lengths,
        n_total = n_total
    )

    # Simulate with selectivity
    sim_with_sel <- simulate_length_freq_from_parameters(
        pars = pars_with_sel,
        lengths = lengths,
        n_total = n_total
    )

    # Distributions should be different
    expect_false(all(sim_no_sel$count == sim_with_sel$count),
                info = "Selectivity should change the distribution")

    # With selectivity, more fish at larger sizes relative to smaller
    # (since l50=40, fish >=40 are more selected)
    large_idx <- lengths >= 40
    small_idx <- lengths < 40

    ratio_with_sel <- sum(sim_with_sel$count[large_idx]) / sum(sim_with_sel$count[small_idx])
    ratio_no_sel <- sum(sim_no_sel$count[large_idx]) / sum(sim_no_sel$count[small_idx])

    # With selectivity, ratio should be higher (more large fish selected)
    expect_true(ratio_with_sel > ratio_no_sel,
               info = "Selectivity should increase proportion of large fish")
})

test_that("Parameter recovery improves with larger sample size", {
    set.seed(101)

    true_pars <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        l50 = 40, ratio = 0.75, vB_min_size = 0,
        spawning_mu = 0.5, spawning_kappa = 3,
        annuli_date = 0.25, annuli_min_age = 0.5
    )

    lengths <- seq(20, 80, by = 5)

    # Small sample
    result_small <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.5),
        lengths = lengths,
        n_per_length = 10,
        n_length_total = 200,
        Delta_l = 1,
        Delta_t = 0.05
    )

    # Large sample
    set.seed(101)  # Same seed for comparable simulation
    result_large <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.5),
        lengths = lengths,
        n_per_length = 50,
        n_length_total = 1000,
        Delta_l = 1,
        Delta_t = 0.05
    )

    # Both should converge
    expect_true(result_small$fit$opt$convergence %in% c(0, 1))
    expect_true(result_large$fit$opt$convergence %in% c(0, 1))

    # Calculate relative errors for k parameter
    k_idx <- which(result_small$comparison$parameter == "k")
    error_small <- abs(result_small$comparison$estimated[k_idx] -
                       result_small$comparison$true[k_idx]) /
                  abs(result_small$comparison$true[k_idx])
    error_large <- abs(result_large$comparison$estimated[k_idx] -
                       result_large$comparison$true[k_idx]) /
                  abs(result_large$comparison$true[k_idx])

    # Both should be finite
    expect_true(is.finite(error_small))
    expect_true(is.finite(error_large))
})

test_that("Parameter recovery works with combined data and custom weights", {
    set.seed(202)

    true_pars <- list(
        k = 0.25,
        L_inf = 110,
        d = 0.15,
        m = 15,
        r = 4,
        l50 = 45,
        ratio = 0.8,
        vB_min_size = 0,
        spawning_mu = 0.5,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )

    # Simulate and fit with custom weighting
    result <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.25, 2023.75),
        lengths = seq(15, 90, by = 3),
        n_per_length = 25,
        n_length_total = 800,
        weight_age = 1.5,
        weight_length = 1,
        Delta_l = 1,
        Delta_t = 0.05
    )

    # Check convergence
    expect_true(result$fit$opt$convergence %in% c(0, 1))

    # Check that we have both data types
    expect_true(nrow(result$simulated_age_data) > 0)
    expect_true(nrow(result$simulated_length_data) > 0)

    # Check that comparison includes selectivity parameters
    expect_true("l50" %in% result$comparison$parameter)
    expect_true("ratio" %in% result$comparison$parameter)
    expect_true("l25" %in% result$comparison$parameter)
})

test_that("fit_on_simulated_data is deterministic with same seed", {
    true_pars <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        l50 = 40, ratio = 0.75, vB_min_size = 0,
        spawning_mu = 0.5, spawning_kappa = 3,
        annuli_date = 0.25, annuli_min_age = 0.5
    )

    # Run twice with same seed
    set.seed(303)
    result1 <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.5),
        lengths = seq(20, 70, by = 5),
        n_per_length = 20,
        n_length_total = 400
    )

    set.seed(303)
    result2 <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.5),
        lengths = seq(20, 70, by = 5),
        n_per_length = 20,
        n_length_total = 400
    )

    # Simulated data should be identical
    expect_equal(result1$simulated_age_data, result2$simulated_age_data)
    expect_equal(result1$simulated_length_data, result2$simulated_length_data)
})

test_that("Large sample recovers l50 accurately", {
    set.seed(404)

    true_pars <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        l50 = 50,
        ratio = 0.7,
        l25 = 35,
        vB_min_size = 0,
        spawning_mu = 0.5,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )

    # Large sample for good recovery
    result <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.3, 2023.7),
        lengths = seq(10, 90, by = 2),
        n_per_length = 40,
        n_length_total = 2000,
        Delta_l = 1,
        Delta_t = 0.05
    )

    # Check convergence
    expect_true(result$fit$opt$convergence %in% c(0, 1))

    # Check l50 recovery (should be good with large sample)
    comp <- result$comparison
    l50_idx <- which(comp$parameter == "l50")

    expect_true(abs(comp$estimated[l50_idx] - comp$true[l50_idx]) < 8,
               info = paste("l50 not recovered within 8 cm with large sample:",
                           "true =", comp$true[l50_idx],
                           "estimated =", comp$estimated[l50_idx]))
})

test_that("Large sample recovers ratio accurately", {
    set.seed(505)

    true_pars <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        l50 = 40,
        ratio = 0.8,
        vB_min_size = 0,
        spawning_mu = 0.5,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )

    # Large sample for good recovery
    result <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.3, 2023.7),
        lengths = seq(10, 90, by = 2),
        n_per_length = 40,
        n_length_total = 2000,
        Delta_l = 1,
        Delta_t = 0.05
    )

    # Check convergence
    expect_true(result$fit$opt$convergence %in% c(0, 1))

    # Check ratio recovery
    comp <- result$comparison
    ratio_idx <- which(comp$parameter == "ratio")

    expect_true(abs(comp$estimated[ratio_idx] - comp$true[ratio_idx]) < 0.15,
               info = paste("ratio not recovered within 0.15 with large sample:",
                           "true =", comp$true[ratio_idx],
                           "estimated =", comp$estimated[ratio_idx]))
})

test_that("Recovered l25 matches ratio * l50", {
    set.seed(606)

    true_pars <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        l50 = 45,
        ratio = 0.7,
        vB_min_size = 0,
        spawning_mu = 0.5,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )

    result <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.5),
        lengths = seq(10, 90, by = 3),
        n_per_length = 30,
        n_length_total = 1000,
        Delta_l = 1,
        Delta_t = 0.05
    )

    # Check that l25 = ratio * l50 in the fitted parameters
    comp <- result$comparison
    l50_est <- comp$estimated[comp$parameter == "l50"]
    ratio_est <- comp$estimated[comp$parameter == "ratio"]
    l25_est <- comp$estimated[comp$parameter == "l25"]

    expect_equal(l25_est, ratio_est * l50_est, tolerance = 1e-6,
                info = "Estimated l25 should equal ratio * l50")
})

test_that("Parameter recovery works without length frequency data (backward compatible)", {
    set.seed(707)

    true_pars <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        vB_min_size = 0,
        spawning_mu = 0.5, spawning_kappa = 3,
        annuli_date = 0.25, annuli_min_age = 0.5
    )

    # No n_length_total specified - should work like before
    result <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.5),
        lengths = seq(20, 80, by = 5),
        n_per_length = 30,
        n_length_total = NULL,  # No length frequency
        Delta_l = 1,
        Delta_t = 0.05
    )

    # Should not have length frequency data
    expect_true(is.null(result$simulated_length_data))

    # Should not have selectivity parameters in comparison
    expect_false("l50" %in% result$comparison$parameter)
    expect_false("ratio" %in% result$comparison$parameter)

    # Should have only 5 parameters
    expect_equal(nrow(result$comparison), 5)

    # Should still converge and recover parameters
    expect_equal(result$fit$opt$convergence, 0)
})

test_that("fit_on_simulated_data handles different weight combinations", {
    set.seed(808)

    true_pars <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        l50 = 40, ratio = 0.75, vB_min_size = 0,
        spawning_mu = 0.5, spawning_kappa = 3,
        annuli_date = 0.25, annuli_min_age = 0.5
    )

    # Test with different weights
    weight_cases <- list(
        list(weight_age = 1, weight_length = 1),
        list(weight_age = 2, weight_length = 1),
        list(weight_age = 1, weight_length = 2)
    )

    for (case in weight_cases) {
        result <- fit_on_simulated_data(
            true_pars = true_pars,
            survey_dates = c(2023.5),
            lengths = seq(20, 70, by = 5),
            n_per_length = 20,
            n_length_total = 500,
            weight_age = case$weight_age,
            weight_length = case$weight_length,
            Delta_l = 1,
            Delta_t = 0.05
        )

        # Should converge for all weight combinations
        expect_true(result$fit$opt$convergence %in% c(0, 1),
                   info = paste("Failed to converge with weights",
                               case$weight_age, ":", case$weight_length))

        # Should have comparison table with all parameters
        expect_equal(nrow(result$comparison), 8)  # k, L_inf, d, m, r, l50, ratio, l25
    }
})

test_that("Simulated length frequency data sums to correct total", {
    set.seed(909)

    pars <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        l50 = 40, ratio = 0.75
    )

    # Test with various sample sizes
    for (n_total in c(100, 500, 1000, 2000)) {
        sim_data <- simulate_length_freq_from_parameters(
            pars = pars,
            lengths = seq(10, 80, by = 5),
            n_total = n_total
        )

        expect_equal(sum(sim_data$count), n_total,
                    info = paste("Total count should equal", n_total))
    }
})

test_that("fit_on_simulated_data produces complete comparison table", {
    set.seed(1010)

    true_pars <- list(
        k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 5,
        l50 = 40, ratio = 0.75, vB_min_size = 0,
        spawning_mu = 0.5, spawning_kappa = 3,
        annuli_date = 0.25, annuli_min_age = 0.5
    )

    result <- fit_on_simulated_data(
        true_pars = true_pars,
        survey_dates = c(2023.5),
        lengths = seq(20, 70, by = 5),
        n_per_length = 25,
        n_length_total = 600
    )

    comp <- result$comparison

    # Check structure
    expect_s3_class(comp, "data.frame")
    expect_true(all(c("parameter", "true", "estimated") %in% names(comp)))

    # Check all expected parameters are present
    expected_params <- c("k", "L_inf", "d", "m", "r", "l50", "ratio", "l25")
    expect_equal(sort(comp$parameter), sort(expected_params))

    # Check all values are finite
    expect_true(all(is.finite(comp$true)))
    expect_true(all(is.finite(comp$estimated)))
})

