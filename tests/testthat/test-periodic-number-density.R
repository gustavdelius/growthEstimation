test_that("periodic density is approximately 1-year periodic", {
    pars <- list(k = 0.2, L_inf = 100, d = 0.1, m = 0.1, r = 0.3,
                 spawning_mu = 0.5, spawning_kappa = 4)

    l_max <- 50
    Delta_l <- 1
    Delta_t <- 0.1
    t_max <- 20

    u_periodic <- getPeriodicNumberDensity(pars, l_max = l_max,
                                           Delta_l = Delta_l,
                                           t_max = t_max,
                                           Delta_t = Delta_t)

    # indices for t and t+1 (t=0..1)
    # grid is t = 0, 0.1, 0.2, ..., 2.0  (N = 21 rows)
    # so t=0 -> 1, t=1 -> 11
    idx_t0 <- 1
    idx_t1 <- 11

    # Compare entire length profiles at t and t+1
    diffs <- abs(u_periodic[idx_t0, ] - u_periodic[idx_t1, ])

    # Allow small numerical tolerance
    expect_lt(max(diffs), 1e-6)
})

