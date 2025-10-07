#' Fit growth parameters by minimizing the negative log likelihood with TMB
#'
#' Optimizes k, L_inf, d, m, r, and optionally l50 and ratio (l25/l50) using nlminb().
#' Spawning parameters, annuli_date, and annuli_min_age are treated as data.
#' Can fit to age-at-length data, length frequency data, or both with weighting.
#'
#' @param pars List of parameters. Must include l50 and ratio if length_freq is provided.
#' @param age_at_length Data frame with columns survey_date (numeric), length, K, count.
#' @param length_freq Optional data frame with columns length and count for length frequency data.
#' @param weight_age Numeric weight for age-at-length likelihood. Default 1.
#' @param weight_length Numeric weight for length frequency likelihood. Default 1.
#' @param Delta_l Numeric size step (cm), default 1.
#' @param Delta_t Numeric time step (years), default 0.05.
#' @param lower Named numeric vector of lower bounds.
#' @param upper Named numeric vector of upper bounds.
#' @return A list with updated `pars`, optimizer result, sdreport, obj, and data
#'   used.
#' @export
fit_tmb_nll <- function(
    pars,
    age_at_length,
    length_freq = NULL,
    weight_age = 1,
    weight_length = 1,
    Delta_l = 1,
    Delta_t = 0.05,
    lower = c(),
    upper = c()
) {
    stopifnot(all(c("survey_date", "length", "K", "count") %in% names(age_at_length)))
    age_at_length <- as.data.frame(age_at_length)
    
    # Ensure r parameter is present (do this early)
    if (is.null(pars$r)) {
        pars$r <- 0  # Default initial value
    }
    
    # Check if length_freq is provided
    use_length_freq <- !is.null(length_freq)
    if (use_length_freq) {
        stopifnot(all(c("length", "count") %in% names(length_freq)))
        length_freq <- as.data.frame(length_freq)
        stopifnot(!is.null(pars$l50), !is.null(pars$ratio))
    }

    # Ensure vB_min_size present in pars ----
    if (is.null(pars$vB_min_size)) {
        pars$vB_min_size <- as.numeric(min(age_at_length$length))
    }

    # Build grids similarly to get_age_log_likelihood()
    l_max <- ceiling(max(age_at_length$length) * 1.1)
    t_max <- max(age_at_length$K) + 2
    N_l <- ceiling(l_max / Delta_l)
    N_t <- ceiling(t_max / Delta_t)
    l_grid <- (1:N_l - 0.5) * Delta_l
    a_grid <- (0:N_t) * Delta_t

    # Map observed length to nearest grid cell index
    length_index <- pmax(1L, pmin(N_l, as.integer(floor(age_at_length$length / Delta_l + 0.5))))

    # Unique age_at_length and indexing
    survey_levels <- sort(unique(age_at_length$survey_date))
    survey_index <- match(age_at_length$survey_date, survey_levels)

    # Observations as vectors, aggregated
    obs_df <- data.frame(
        s = survey_index,
        j = length_index,
        K = as.integer(age_at_length$K),
        count = as.numeric(age_at_length$count)
    )
    obs_df <- stats::aggregate(count ~ s + j + K, data = obs_df, sum)
    obs_survey_index <- as.integer(obs_df$s)
    obs_length_index <- as.integer(obs_df$j)
    obs_K <- as.integer(obs_df$K)
    obs_count <- as.numeric(obs_df$count)

    # Prepare length frequency data
    if (use_length_freq) {
        # Map length_freq lengths to grid indices
        length_freq_index <- pmax(1L, pmin(N_l, as.integer(ceiling(length_freq$length / Delta_l))))
        length_freq_count <- as.numeric(length_freq$count)
    } else {
        # Provide dummy data if not using length frequency
        length_freq_index <- integer(0)
        length_freq_count <- numeric(0)
    }

    # Prepare TMB data and parameters
    tmb_data <- list(
        obs_survey_index = obs_survey_index,
        obs_length_index = obs_length_index,
        obs_K = obs_K,
        obs_count = obs_count,
        length_freq_index = length_freq_index,
        length_freq_count = length_freq_count,
        survey_dates = as.numeric(survey_levels),
        l_grid = as.numeric(l_grid),
        a_grid = as.numeric(a_grid),
        spawning_mu = as.numeric(pars$spawning_mu),
        spawning_kappa = as.numeric(pars$spawning_kappa),
        annuli_date = as.numeric(pars$annuli_date),
        annuli_min_age = as.numeric(pars$annuli_min_age),
        Delta_l = as.numeric(Delta_l),
        Delta_t = as.numeric(Delta_t),
        log_eps = log(1e-9),
        vB_min_size = as.numeric(pars$vB_min_size),
        weight_age = as.numeric(weight_age),
        weight_length = as.numeric(weight_length)
    )

    # Define parameter names and extract based on whether we're using length frequency
    if (use_length_freq) {
        parameter_names <- c("k", "L_inf", "d", "m", "r", "l50", "ratio")
        tmb_parameters <- list(
            k = pars$k,
            L_inf = pars$L_inf,
            d = pars$d,
            m = pars$m,
            r = pars$r,
            l50 = pars$l50,
            ratio = pars$ratio
        )
    } else {
        parameter_names <- c("k", "L_inf", "d", "m", "r")
        tmb_parameters <- list(
            k = pars$k,
            L_inf = pars$L_inf,
            d = pars$d,
            m = pars$m,
            r = pars$r,
            l50 = 40.0,   # Dummy value (will be fixed via map)
            ratio = 0.75  # Dummy value (will be fixed via map)
        )
    }

    # Set up bounds
    if (use_length_freq) {
        lower_limit = c(k = 1e-6, L_inf = 1e-3, d = 1e-6, m = 1e-6, r = 1e-6, 
                       l50 = 1e-3, ratio = 1e-6)
        upper_limit = c(k = Inf, L_inf = Inf, d = Inf, m = Inf, r = Inf, 
                       l50 = Inf, ratio = 0.999)  # ratio < 1 to enforce l25 < l50
    } else {
        lower_limit = c(k = 1e-6, L_inf = 1e-3, d = 1e-6, m = 1e-6, r = 1e-6)
        upper_limit = c(k = Inf, L_inf = Inf, d = Inf, m = Inf, r = Inf)
    }
    
    if (!all(names(lower) %in% parameter_names)) {
        bad_names <- setdiff(names(lower), parameter_names)
        stop("You cannot specify a lower limit on: ", paste(bad_names, collapse = ", "))
    }
    lower_limit[names(lower)] <- lower[names(lower)]

    if (!all(names(upper) %in% parameter_names)) {
        bad_names <- setdiff(names(upper), parameter_names)
        stop("You cannot specify an upper limit on: ", paste(bad_names, collapse = ", "))
    }
    upper_limit[names(upper)] <- upper[names(upper)]

    # Create map to fix parameters if not using length frequency
    if (use_length_freq) {
        obj <- TMB::MakeADFun(
            data = tmb_data,
            parameters = tmb_parameters,
            DLL = "growthEstimation",
            silent = TRUE
        )
    } else {
        # Fix l50 and ratio when not using length frequency
        map_list <- list(l50 = factor(NA), ratio = factor(NA))
        obj <- TMB::MakeADFun(
            data = tmb_data,
            parameters = tmb_parameters,
            map = map_list,
            DLL = "growthEstimation",
            silent = TRUE
        )
    }

    opt <- nlminb(
        start = obj$par,
        objective = obj$fn,
        gradient = obj$gr,
        lower = as.numeric(lower_limit[parameter_names]),
        upper = as.numeric(upper_limit[parameter_names])
    )

    # Update parameters in `pars`
    par <- opt$par
    pars[names(par)] <- par[names(par)]
    
    # If using length frequency, convert ratio back to l25
    if (use_length_freq) {
        pars$l25 <- pars$ratio * pars$l50
    }

    sdr <- try(TMB::sdreport(obj), silent = TRUE)

    list(
        pars = pars,
        opt = opt,
        sdreport = sdr,
        obj = obj,
        data = tmb_data
    )
}



