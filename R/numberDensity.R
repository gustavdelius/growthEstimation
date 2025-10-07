#' von Mises Probability Density Function
#'
#' A lightweight implementation used to model circular seasonality
#' (e.g., spawning day within a year).
#' @param x Angle in radians.
#' @param mu Mean direction in radians.
#' @param kappa Concentration parameter (higher means more concentrated).
#' @return The probability density evaluated at `x`.
#' @export
#' @examples
#' # Density at angle pi when centered at pi with moderate concentration
#' von_mises_pdf(pi, mu = pi, kappa = 2)
von_mises_pdf <- function(x, mu, kappa) {
    denominator <- 2 * pi * besselI(kappa, 0)
    numerator <- exp(kappa * cos(x - mu))
    return(numerator / denominator)
}

#' Spawning Density Function
#'
#' Calculates relative spawning intensity S(d) for numeric dates using a
#' von Mises density on the unit circle of the year.
#' @param numeric_dates A vector of numeric dates (e.g., 2023.45). Only the
#'   fractional part is used to determine day-of-year.
#' @param mu The mean spawning date as a fraction of a year in \[0, 1). For
#'   example, 0.5 is mid-year.
#' @param kappa The concentration parameter for the spawning distribution.
#' @return A numeric vector of relative spawning intensities with the same
#'   length as `numeric_dates`.
#' @export
#' @examples
#' # Peak at mid-year with moderate spread
#' get_spawning_density(numeric_dates = c(2020.45, 2020.50, 2020.55),
#'                  mu = 0.5, kappa = 4)
get_spawning_density <- function(numeric_dates, mu, kappa) {
    day_fraction <- numeric_dates %% 1
    day_rad <- day_fraction * 2 * pi
    mu_rad <- mu * 2 * pi
    density <- von_mises_pdf(day_rad, mu = mu_rad, kappa = kappa)
    return(density)
}

#' Get number density as a function of length and time
#'
#' Convolves the Green's function obtained with `get_greens_function()` with the
#' spawning density from `get_spawning_density()` to produce the number density
#' as a function of length and time.
#' @inheritParams get_greens_function
#' @param single_cohort Logical. If TRUE (default) include only individuals born
#'   within the first year. If FALSE, includes all individuals born between time
#'   0 and time `t_max`.
#' @return A matrix holding the number density u(t,l). Rows correspond to time
#'   and columns to length.
#' @export
get_number_density <- function(pars, l_max, Delta_l = 1,
                               t_max = 10, Delta_t = 0.05,
                               single_cohort = TRUE) {
    G <- get_greens_function(pars, l_max = l_max, Delta_l = Delta_l,
                             t_max = t_max, Delta_t = Delta_t)
    # G has rows = time (age) 0..N_t and cols = size classes

    # Build age/time grid
    N_t <- nrow(G) - 1L
    a_grid <- (0:N_t) * Delta_t

    # Output matrix u with same shape as G (rows time, cols size)
    u <- matrix(0, nrow = N_t + 1L, ncol = ncol(G))

    # Spawning parameters
    mu <- pars$spawning_mu
    kappa <- pars$spawning_kappa

    # For each time t_n, convolve S(t_n - a) with G(a, l) over a in [0, t_n]
    for (n in 0:N_t) {
        # ages considered up to current time
        idx_all <- 0:n
        if (single_cohort) {
            # keep only ages in (t-1, t] (restrict to >=0 automatically)
            t_n <- a_grid[n + 1L]
            lower <- max(0, t_n - 1)
            ages <- a_grid[idx_all + 1L]
            keep <- (ages > lower) & (ages <= t_n)
            idx <- idx_all[keep]
        } else {
            idx <- idx_all
        }
        if (length(idx) == 0L) {
            next
        }
        # birth dates as numeric years (fractional part encodes day-of-year)
        birth_dates <- a_grid[n + 1L] - a_grid[idx + 1L]
        w <- get_spawning_density(birth_dates, mu = mu, kappa = kappa)
        # Weighted sum over ages (rows) to get u(t_n, l)
        # Using crossprod to get 1 x N_l row
        u[n + 1L, ] <- as.numeric(crossprod(w, G[idx + 1L, , drop = FALSE])) * Delta_t
    }

    return(u)
}

#' Get periodic number density as a function of length and time
#'
#' Calculates the population density over time under the assumption that births
#' have been happening for all negative years in the past, resulting in a
#' perfectly periodic density with a period of 1 year.
#'
#' This function extends the convolution approach used in `get_number_density()`
#' to include contributions from all past years by exploiting the periodicity of
#' the spawning density function.
#'
#' @inheritParams get_number_density
#' @return A matrix holding the periodic number density u(t,l). Rows correspond
#'   to time and columns to length. The density is periodic with period 1 year.
#' @export
#' @examples
#' # Get periodic density for a simple parameter set
#' pars <- list(k = 0.2, L_inf = 100, d = 0.1, m = 0.1,
#'              spawning_mu = 0.5, spawning_kappa = 4)
#' u_periodic <- get_periodic_number_density(pars, l_max = 50, t_max = 2)
get_periodic_number_density <- function(pars, l_max, Delta_l = 1,
                                        t_max = 1, Delta_t = 0.05) {
    # For periodic solution, we need a Green's function that covers enough time
    # so that the oldest cohorts have essentially died out.
    # Use a larger t_max for Green's function to ensure we capture all
    # contributions. At least 3x the simulation time, minimum 20 years.
    greens_t_max <- max(t_max * 3, 20)

    G <- get_greens_function(pars, l_max = l_max, Delta_l = Delta_l,
                   t_max = greens_t_max, Delta_t = Delta_t)
    # G has rows = time (age) 0..N_t_greens and cols = size classes

    max_attempts <- 3
    attempt <- 1
    while (!all(G[nrow(G), ] < 1e-4)) {
        if (attempt >= max_attempts) {
            stop("Cohorts live for too long (density at max age > 1e-4 after ",
                 greens_t_max,
                 " years). Please increase the mortality rate (pars$m).")
        }
        greens_t_max <- greens_t_max * 2
        message("I am running the cohorts for ", greens_t_max, " years.")
        G <- get_greens_function(pars,
            l_max = l_max, Delta_l = Delta_l,
            t_max = greens_t_max, Delta_t = Delta_t
        )
        attempt <- attempt + 1
    }

    # Build age/time grid for Green's function
    N_t_greens <- nrow(G) - 1L
    a_grid <- (0:N_t_greens) * Delta_t

    # Build time grid for output (only up to requested t_max)
    N_t_output <- ceiling(t_max / Delta_t)
    t_grid <- (0:N_t_output) * Delta_t

    # Output matrix u (rows = output time, cols = size)
    u <- matrix(0, nrow = N_t_output + 1L, ncol = ncol(G))

    # Spawning parameters
    mu <- pars$spawning_mu
    kappa <- pars$spawning_kappa

    # For periodic solution, we need to consider contributions from all past
    # years The key insight is that
    # get_spawning_density(t) = get_spawning_density(t + 1) for any t.
    # So we can sum over all integer years in the past.

    # For each time t_n, we need to sum over all ages a, including those
    # that would correspond to births in negative years.
    for (n in 0:N_t_output) {
        t_n <- t_grid[n + 1L]

        # For each possible age a (including ages > t_n for periodic solution)
        # We need to consider all ages up to the maximum age in the Green's
        # function
        for (a_idx in 0:N_t_greens) {
            a <- a_grid[a_idx + 1L]

            # Birth date for this age cohort
            birth_date <- t_n - a

            # For periodic solution, we need to sum over all integer years
            # The spawning density is periodic with period 1, so:
            # S(birth_date) = S(birth_date + k) for any integer k

            # We can calculate the contribution from this age cohort
            # by evaluating the spawning density at the fractional part of
            # birth_date
            birth_date_fractional <- birth_date %% 1
            spawning_intensity <- get_spawning_density(birth_date_fractional, mu = mu, kappa = kappa)

            # Add contribution from this age cohort
            u[n + 1L, ] <- u[n + 1L, ] +
                spawning_intensity * G[a_idx + 1L, ] * Delta_t
        }
    }

    return(u)
}

#' Get log likelihood of length frequency observations
#'
#' Calculates the multinomial log likelihood of observed length frequencies
#' given the steady state density predicted by the model, accounting for
#' size selectivity.
#'
#' @param pars A list containing the model parameters to pass to
#'   `solve_pde_steady_state()`: k, L_inf, d, m, r (and optionally vB_min_size),
#'   as well as selectivity parameters l50 and l25.
#' @param length_freq A data frame with columns `length` and `count`, containing
#'   observed length frequencies.
#' @param Delta_l Width of size bins (cm). Default is 1.
#' @param l_max The maximum length to consider. If NULL (default), it is set to
#'   110% of the maximum observed length.
#'
#' @return A scalar value representing the log likelihood (negative of the
#'   multinomial negative log likelihood).
#' @export
#' @examples
#' # Example using the Cod Celtic Sea data
#' data(Cod_CS_pars)
#' data(Cod_CS_length_freq)
#' # Add selectivity parameters
#' pars <- Cod_CS_pars
#' pars$l50 <- 40
#' pars$l25 <- 30
#' get_length_log_likelihood(pars, Cod_CS_length_freq)
get_length_log_likelihood <- function(pars, length_freq, Delta_l = 1, l_max = NULL) {
    
    # Determine l_max if not provided
    if (is.null(l_max)) {
        l_max <- ceiling(max(length_freq$length) * 1.1)
    }
    
    # Get steady state density
    u_steady <- solve_pde_steady_state(pars, Delta_l = Delta_l, l_max = l_max)
    
    # Create length grid (cell centers)
    N_l <- length(u_steady)
    l_grid <- (1:N_l - 0.5) * Delta_l
    
    # Apply size selectivity if parameters are provided
    if (!is.null(pars$l50) && !is.null(pars$l25)) {
        # Calculate slope from l25 and l50
        # At l50: selectivity = 0.5, at l25: selectivity = 0.25
        # Using logistic function: S(l) = 1 / (1 + exp(-slope * (l - l50)))
        # Solving: slope = log(3) / (l50 - l25)
        slope <- log(3) / (pars$l50 - pars$l25)
        
        # Calculate selectivity for each length in the grid
        selectivity <- 1 / (1 + exp(-slope * (l_grid - pars$l50)))
        
        # Multiply density by selectivity
        u_steady <- u_steady * selectivity
    }
    
    # Normalize density to get probabilities
    total_density <- sum(u_steady)
    if (total_density > 0) {
        prob_model <- u_steady / total_density
    } else {
        # If density is zero everywhere, use uniform distribution
        prob_model <- rep(1 / N_l, N_l)
    }
    
    # Match observed lengths to grid
    # For each observed length, find the corresponding grid cell
    length_indices <- pmax(1, pmin(N_l, ceiling(length_freq$length / Delta_l)))
    
    # Get probabilities for observed lengths
    prob_obs <- prob_model[length_indices]
    
    # Calculate multinomial log likelihood
    # Add small epsilon to avoid log(0)
    epsilon <- 1e-9
    log_lik <- sum(length_freq$count * log(prob_obs + epsilon))
    
    return(log_lik)
}
