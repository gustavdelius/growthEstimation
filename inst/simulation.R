# Here we simulate age-at-size data and see how well we can
# estimate the parameters from the simulated data.

true_pars <- list(
    k = 0.2, L_inf = 90, d = 0.5, m = 4, r = 5,
    l50 = 20, l25 = 15,
    spawning_mu = 0.5, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
)
start_pars <- list(
    k = 0.1, L_inf = 150, d = 0.1, m = 10, r = 1,
    l50 = 10, l25 = 9,
    spawning_mu = 0.5, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
)

survey_dates <- c(2023.25, 2023.5, 2023.75)
lengths <- seq(10, 100, by = 1)
n_per_length = 20
n_total <- 5000

set.seed(40)
res <- fit_on_simulated_data(
    true_pars = true_pars,
    start_pars = start_pars,
    survey_dates = survey_dates,
    lengths = lengths,
    n_per_length = n_per_length,
    n_length_total = n_total,
    weight_length = 6,
    lower = c()
)
res$fit$opt$convergence
res$comparison

plot_age_likelihood(res$fit$pars, res$simulated_age_data)
plot_age_likelihood(true_pars, res$simulated_age_data)
plot_length(res$fit$pars, res$simulated_length_data)
plot_length(true_pars, res$simulated_length_data)

u <- get_number_density(res$fit$pars, l_max = 100)
plot_density_2d(u)
plot_density_3d(u, l_min = 10)


