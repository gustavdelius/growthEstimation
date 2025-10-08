#' Get log likelihood of age observations
#'
#' @param pars A list containing the model parameters: k, L_inf, d, m, r.
#' @param age_at_length A data frame with survey age-at-length observations with
#'   columns `survey_date`, `length`, `K`, and `count`.
#' @param Delta_l Width of size bins (cm). Default is 1.
#' @param Delta_t Time step for the model simulation (years). Default is 0.05.
#'
#' @return A data frame with, for each observed length-K bin in each survey, the
#'   observed count, expected count under the model, model probability, sample
#'   size, negative log-likelihood contribution, and signed negative
#'   log-likelihood contribution.
#' @export
get_age_log_likelihood <- function(pars, age_at_length, Delta_l = 1, Delta_t = 0.05) {

    # Set up grids
    l_max <- ceiling(max(age_at_length$length) * 1.1) # A bit larger than maximum observed size
    t_max <- max(age_at_length$K) + 2 # TODO: set t_max only as high as needed
    N_l <- ceiling(l_max / Delta_l) # Number of time steps
    N_t <- ceiling(t_max / Delta_t) # Number of time steps
    l_grid <- (1:N_l) * Delta_l # Size grid (bin start values)
    a_grid <- (0:N_t) * Delta_t   # Age grid (including age = 0)

    # Ensure vB_min_size in pars once for downstream calls
    if (is.null(pars$vB_min_size)) {
        pars$vB_min_size <- as.numeric(min(age_at_length$length))
    }
    G <- get_greens_function(pars, l_max = l_max, Delta_l = Delta_l,
                   t_max = t_max, Delta_t = Delta_t)

    age_likelihood(
        age_at_length, G = G, a = a_grid, l = l_grid,
        mu = pars$spawning_mu, kappa = pars$spawning_kappa,
        annuli_date = pars$annuli_date, annuli_min_age = pars$annuli_min_age)
}


"Simulate survey data and parameter recovery helpers"

#' Simulate age-at-length survey observations from model parameters
#'
#' Generates multinomial samples of annuli counts K for specified age_at_length and
#' length distributions, using the model predictions implied by `pars`.
#'
#' @param pars List with at least: k, L_inf, d, m, r, spawning_mu,
#'   spawning_kappa, annuli_date, annuli_min_age. Optionally `vB_min_size`.
#' @param survey_dates Numeric vector of survey dates (e.g. c(2023.25, 2023.75)).
#' @param lengths Numeric vector of length-class centers to observe.
#' @param n_per_length Integer or vector: number of fish sampled per length per
#'   survey. If a single integer, it is recycled for all lengths and age_at_length.
#' @param Delta_l Numeric size step (cm). Default 1.
#' @param Delta_t Numeric time step (years). Default 0.05.
#' @return Data frame with columns `survey_date`, `length`, `K`, `count`.
#' @export
simulate_age_at_length_from_parameters <- function(
    pars,
    survey_dates,
    lengths,
    n_per_length,
    Delta_l = 1,
    Delta_t = 0.05
) {
    stopifnot(length(survey_dates) >= 1)
    stopifnot(length(lengths) >= 1)

    # Ensure vB_min_size present
    if (is.null(pars$vB_min_size)) {
        pars$vB_min_size <- as.numeric(min(lengths))
    }

    # Build grids consistent with model likelihood
    l_max <- ceiling(max(lengths) * 1.1)
    # Reasonable upper age to cover plausible K
    t_max <- 10
    N_l <- ceiling(l_max / Delta_l)
    l_grid <- (1:N_l) * Delta_l
    a_grid <- (0:ceiling(t_max / Delta_t)) * Delta_t

    # Greens function once (independent of survey)
    G <- get_greens_function(pars, l_max = l_max, Delta_l = Delta_l,
                   t_max = t_max, Delta_t = Delta_t)

    # Map desired observation lengths to nearest model grid rows for P(K|l)
    # We'll obtain P_model on model grid then subset by nearest length classes
    len_indices <- pmax(1L, pmin(N_l, as.integer(floor(lengths / Delta_l + 0.5))))
    obs_lengths <- (len_indices) * Delta_l

    # Handle n_per_length vectorization
    if (length(n_per_length) == 1L) {
        n_len <- rep(as.integer(n_per_length), length(lengths))
    } else {
        stopifnot(length(n_per_length) == length(lengths))
        n_len <- as.integer(n_per_length)
    }

    results <- list()
    res_i <- 1L
    for (sd in survey_dates) {
        # Get P(K|l) for this survey date
        P_model <- annuli_predictions_for_date(
            survey_date = sd,
            G = G,
            a = a_grid,
            l = l_grid,
            mu = pars$spawning_mu,
            kappa = pars$spawning_kappa,
            annuli_date = pars$annuli_date,
            annuli_min_age = pars$annuli_min_age
        )

        # Subset to observed lengths
        P_obs <- P_model[as.character(obs_lengths), , drop = FALSE]

        # Simulate K counts per length via multinomial
        k_vals <- as.numeric(colnames(P_obs))
        for (ii in seq_along(obs_lengths)) {
            probs <- P_obs[ii, ]
            # Guard against zero-prob rows
            if (!is.finite(sum(probs)) || sum(probs) == 0) {
                next
            }
            counts <- as.vector(rmultinom(1, size = n_len[ii], prob = probs))
            if (sum(counts) > 0) {
                results[[res_i]] <- data.frame(
                    survey_date = sd,
                    length = lengths[ii],
                    K = k_vals,
                    count = counts
                )
                res_i <- res_i + 1L
            }
        }
    }

    if (length(results) == 0L) {
        return(data.frame(survey_date = numeric(0), length = numeric(0), K = integer(0), count = numeric(0)))
    }
    out <- do.call(rbind, results)
    rownames(out) <- NULL
    return(out)
}

#' Simulate length frequency data from steady state with selectivity
#'
#' Generates multinomial samples of length frequencies from the steady-state
#' distribution with size selectivity.
#'
#' @param pars List with parameters including k, L_inf, d, m, r, and optionally
#'   l50, l25 or ratio for selectivity.
#' @param lengths Numeric vector of length-class centers to observe.
#' @param n_total Total number of fish to sample.
#' @param Delta_l Numeric size step (cm). Default 1.
#' @param l_max Maximum length for steady state calculation. If NULL, uses 110%
#'   of max(lengths).
#' @return Data frame with columns `length` and `count`.
#' @export
simulate_length_freq_from_parameters <- function(
    pars,
    lengths,
    n_total,
    Delta_l = 1,
    l_max = NULL
) {
    # Calculate l_max if not provided
    if (is.null(l_max)) {
        l_max <- ceiling(max(lengths) * 1.1)
    }

    # Get steady state density
    u_steady <- solve_pde_steady_state(pars, Delta_l = Delta_l, l_max = l_max)

    # Create length grid
    N_l <- length(u_steady)
    l_grid <- (1:N_l - 0.5) * Delta_l

    # Apply selectivity if parameters provided
    if (!is.null(pars$l50)) {
        if (is.null(pars$l25) && !is.null(pars$ratio)) {
            pars$l25 <- pars$ratio * pars$l50
        }

        if (!is.null(pars$l25)) {
            slope <- log(3) / (pars$l50 - pars$l25)
            selectivity <- 1 / (1 + exp(-slope * (l_grid - pars$l50)))
            u_steady <- u_steady * selectivity
        }
    }

    # Normalize to probabilities
    prob_model <- u_steady / sum(u_steady)

    # Map observed lengths to grid indices
    length_indices <- pmax(1, pmin(N_l, ceiling(lengths / Delta_l)))

    # Get probabilities for observed lengths
    prob_obs <- prob_model[length_indices]
    prob_obs <- prob_obs / sum(prob_obs)  # Renormalize to observed lengths

    # Sample from multinomial
    counts <- as.vector(rmultinom(1, size = n_total, prob = prob_obs))

    # Return data frame
    data.frame(
        length = lengths,
        count = counts
    )
}

#' Fit on simulated survey data and compare to true parameters
#'
#' Runs `fit_tmb_nll()` on data from `simulate_age_at_length_from_parameters()` and
#' optionally `simulate_length_freq_from_parameters()`, returning a compact
#' comparison table with true vs estimated parameters.
#'
#' @param true_pars List of true parameters to simulate from.
#' @param survey_dates Numeric vector of survey dates.
#' @param lengths Numeric vector of length-class centers.
#' @param n_per_length Integer or vector sample sizes per length per survey (for age data).
#' @param n_length_total Optional total sample size for length frequency data. If NULL,
#'   no length frequency data is simulated.
#' @param weight_age Weight for age-at-length likelihood. Default 1.
#' @param weight_length Weight for length frequency likelihood. Default 1.
#' @param start_pars Optional list of starting values for optimization; defaults
#'   to `true_pars` with small jitter.
#' @param Delta_l Numeric size step.
#' @param Delta_t Numeric time step.
#' @param lower,upper Optional bound vectors for `fit_tmb_nll()`.
#' @return A list with `simulated_age_data`, `simulated_length_data`, `fit`, and
#'   `comparison` data.frame.
#' @export
fit_on_simulated_data <- function(
    true_pars,
    survey_dates,
    lengths,
    n_per_length,
    n_length_total = NULL,
    weight_age = 1,
    weight_length = 1,
    start_pars = NULL,
    Delta_l = 1,
    Delta_t = 0.05,
    lower = c(),
    upper = c()
) {
    # Simulate age-at-length data
    sim_age_df <- simulate_age_at_length_from_parameters(
        pars = true_pars,
        survey_dates = survey_dates,
        lengths = lengths,
        n_per_length = n_per_length,
        Delta_l = Delta_l,
        Delta_t = Delta_t
    )

    if (nrow(sim_age_df) == 0L) {
        stop("Age simulation produced no observations; check parameters and sampling design.")
    }

    # Simulate length frequency data if requested
    sim_length_df <- NULL
    use_length_freq <- !is.null(n_length_total) && n_length_total > 0

    if (use_length_freq) {
        sim_length_df <- simulate_length_freq_from_parameters(
            pars = true_pars,
            lengths = lengths,
            n_total = n_length_total,
            Delta_l = Delta_l
        )
    }

    # Starting values
    if (is.null(start_pars)) {
        jitter <- function(x) as.numeric(x) * exp(rnorm(1, sd = 0.1))
        start_pars <- true_pars
        param_names <- c("k", "L_inf", "d", "m", "r")
        if (use_length_freq) {
            param_names <- c(param_names, "l50", "ratio")
        }
        for (nm in param_names) {
            if (!is.null(start_pars[[nm]])) start_pars[[nm]] <- jitter(start_pars[[nm]])
        }
    }

    # Fit
    fit <- fit_tmb_nll(
        pars = start_pars,
        age_at_length = sim_age_df,
        length_freq = sim_length_df,
        weight_age = weight_age,
        weight_length = weight_length,
        Delta_l = Delta_l,
        Delta_t = Delta_t,
        lower = lower,
        upper = upper
    )

    # Comparison table
    est_pars <- fit$pars
    keys <- c("k", "L_inf", "d", "m", "r")
    if (use_length_freq) {
        keys <- c(keys, "l50", "ratio", "l25")
    }

    comp <- data.frame(
        parameter = keys,
        true = vapply(keys, function(k) {
            val <- true_pars[[k]]
            if (is.null(val)) NA_real_ else as.numeric(val)
        }, numeric(1)),
        estimated = vapply(keys, function(k) {
            val <- est_pars[[k]]
            if (is.null(val)) NA_real_ else as.numeric(val)
        }, numeric(1)),
        start = vapply(keys, function(k) {
            val <- start_pars[[k]]
            if (is.null(val)) NA_real_ else as.numeric(val)
        }, numeric(1))
    )

    list(
        simulated_age_data = sim_age_df,
        simulated_length_data = sim_length_df,
        fit = fit,
        comparison = comp
    )
}
