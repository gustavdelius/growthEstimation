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
