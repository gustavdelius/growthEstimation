# Here we simulate age-at-size data and see how well we can
# estimate the parameters from the simulated data.

true_pars <- list(
    k = 0.2, L_inf = 90, d = 0.5, m = 4, r = 5,
    spawning_mu = 0.5, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5,
    vB_min_size = 0
)
start_pars <- list(
    k = 0.1, L_inf = 150, d = 0.1, m = 10, r = 1,
    spawning_mu = 0.5, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5,
    vB_min_size = 0
)

u <- getNumberDensity(true_pars, l_max = 100)
plotDensity2D(u)
plotDensity3D(u, l_min = 10)

set.seed(123)
survey_dates <- c(2023.25, 2023.5, 2023.75)
lengths <- seq(10, 100, by = 1)
n_per_length = 50

sim_df <- simulate_surveys_from_parameters(
    pars = true_pars,
    survey_dates = survey_dates,
    lengths = lengths,
    n_per_length = n_per_length
)
plotAgeLikelihood(true_pars, sim_df)

set.seed(123)
res <- fit_on_simulated_data(
    true_pars = true_pars,
    start_pars = start_pars,
    survey_dates = survey_dates,
    lengths = lengths,
    n_per_length = n_per_length
)
res$comparison

u <- getNumberDensity(res$fit$pars, l_max = 100)
plotDensity2D(u)
plotDensity3D(u, l_min = 10)

plotAgeLikelihood(res$fit$pars, sim_df)
