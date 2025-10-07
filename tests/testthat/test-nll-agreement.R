test_that("TMB nll matches NLL computed from get_age_log_likelihood()", {
    age_at_length <- Cod_CS_age_at_length

    # Common discretisation to keep both paths consistent
    Delta_l <- 1
    Delta_t <- 0.05

    # Choose fixed parameters (not optimizing) to make the comparison deterministic
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

    # 1) NLL via get_age_log_likelihood(): sum of NegLogLik contributions
    ll_df <- get_age_log_likelihood(pars, age_at_length, Delta_l = Delta_l, Delta_t = Delta_t)
    nll_r <- sum(ll_df$TotalNegLogLik, na.rm = TRUE)

    # 2) NLL via TMB objective function in nll.cpp using identical grids and data
    tmb <- fit_tmb_nll(
        pars = pars,
        age_at_length = age_at_length,
        Delta_l = Delta_l,
        Delta_t = Delta_t
    )

    # Evaluate the TMB objective at the specified parameters (already set as start)
    nll_tmb <- as.numeric(tmb$obj$fn(tmb$obj$par))

    # Compare within a small tolerance because of numerical details and epsilon terms
    expect_true(is.finite(nll_r) && is.finite(nll_tmb))
    expect_equal(nll_tmb, nll_r, tolerance = 1e-6)
})

test_that("TMB nll matches NLL computed from get_length_log_likelihood()", {
    length_freq <- Cod_CS_length_freq
    age_at_length <- Cod_CS_age_at_length
    
    # Common discretisation
    Delta_l <- 1
    Delta_t <- 0.05
    
    # Fixed parameters with selectivity
    pars <- list(
        k = 0.3,
        L_inf = 100,
        d = 0.2,
        m = 20,
        r = 0.5,
        l50 = 40,
        ratio = 0.75,  # l25/l50
        vB_min_size = 0,  # Important: must match solve_pde_steady_state default
        spawning_mu = 0.4,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )
    
    # Calculate l_max consistently (same as fit_tmb_nll does)
    l_max <- ceiling(max(age_at_length$length, length_freq$length) * 1.1)
    
    # 1) NLL via get_length_log_likelihood()
    # Note: This returns log likelihood (positive contribution to likelihood)
    # so we negate it to get NLL
    log_lik_r <- get_length_log_likelihood(pars, length_freq, Delta_l = Delta_l, l_max = l_max)
    nll_length_r <- -log_lik_r
    
    # 2) NLL via TMB with length frequency only (weight_age = 0, weight_length = 1)
    tmb <- fit_tmb_nll(
        pars = pars,
        age_at_length = age_at_length,
        length_freq = length_freq,
        weight_age = 0,
        weight_length = 1,
        Delta_l = Delta_l,
        Delta_t = Delta_t
    )
    
    # Re-evaluate at starting parameters (before optimization)
    par_start <- c(k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5, l50 = 40, ratio = 0.75)
    tmb$obj$fn(par_start)  # Evaluate to update internal state
    rep <- tmb$obj$report()
    nll_tmb <- rep$nll_length  # Use reported value at starting parameters
    
    # Compare
    expect_true(is.finite(nll_length_r) && is.finite(nll_tmb))
    expect_equal(nll_tmb, nll_length_r, tolerance = 1e-2)  # Slightly larger tolerance for numerical differences
})

test_that("TMB nll matches combined NLL with weighted likelihoods", {
    age_at_length <- Cod_CS_age_at_length
    length_freq <- Cod_CS_length_freq
    
    # Common discretisation
    Delta_l <- 1
    Delta_t <- 0.05
    
    # Fixed parameters
    pars <- list(
        k = 0.3,
        L_inf = 100,
        d = 0.2,
        m = 20,
        r = 0.5,
        l50 = 40,
        ratio = 0.75,
        vB_min_size = 0,  # Important: must match solve_pde_steady_state default
        spawning_mu = 0.4,
        spawning_kappa = 3,
        annuli_date = 0.25,
        annuli_min_age = 0.5
    )
    
    # Weights
    weight_age <- 2.5
    weight_length <- 1.0
    
    # Calculate l_max consistently
    l_max <- ceiling(max(age_at_length$length, length_freq$length) * 1.1)
    
    # 1) NLL from R functions
    ll_age_df <- get_age_log_likelihood(pars, age_at_length, Delta_l = Delta_l, Delta_t = Delta_t)
    nll_age_r <- sum(ll_age_df$TotalNegLogLik, na.rm = TRUE)
    
    log_lik_length_r <- get_length_log_likelihood(pars, length_freq, Delta_l = Delta_l, l_max = l_max)
    nll_length_r <- -log_lik_length_r
    
    nll_combined_r <- weight_age * nll_age_r + weight_length * nll_length_r
    
    # 2) NLL via TMB with weighted likelihoods
    tmb <- fit_tmb_nll(
        pars = pars,
        age_at_length = age_at_length,
        length_freq = length_freq,
        weight_age = weight_age,
        weight_length = weight_length,
        Delta_l = Delta_l,
        Delta_t = Delta_t
    )
    
    # Re-evaluate at starting parameters (before optimization)
    par_start <- c(k = 0.3, L_inf = 100, d = 0.2, m = 20, r = 0.5, l50 = 40, ratio = 0.75)
    tmb$obj$fn(par_start)  # Evaluate to update internal state
    rep <- tmb$obj$report()
    nll_tmb <- rep$nll  # Use reported combined nll at starting parameters
    
    # Compare
    expect_true(is.finite(nll_combined_r) && is.finite(nll_tmb))
    expect_equal(nll_tmb, nll_combined_r, tolerance = 1e-2)  # Slightly larger tolerance for numerical differences
})
