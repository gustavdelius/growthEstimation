#' Fit growth parameters by minimizing the negative log likelihood with TMB
#'
#' Optimizes k, L_inf, d, and m using nlminb().
#' Spawning parameters, annuli_date, and annuli_min_age are treated as data.
#'
#' @param pars List of parameters
#' @param surveys Data frame with columns survey_date (numeric), Length, K, count.
#' @param Delta_l Numeric size step (cm), default 1.
#' @param Delta_t Numeric time step (years), default 0.05.
#' @param lower Named numeric vector of lower bounds.
#' @param upper Named numeric vector of upper bounds.
#' @return A list with updated `pars`, optimizer result, sdreport, obj, and data
#'   used.
#' @export
fit_tmb_nll <- function(
    pars,
    surveys,
    Delta_l = 1,
    Delta_t = 0.05,
    lower = c(),
    upper = c()
) {
    stopifnot(all(c("survey_date", "Length", "K", "count") %in% names(surveys)))
    surveys <- as.data.frame(surveys)

    # Build grids similarly to getLogLik()
    l_max <- ceiling(max(surveys$Length) * 1.1)
    t_max <- max(surveys$K) + 2
    N_l <- ceiling(l_max / Delta_l)
    N_t <- ceiling(t_max / Delta_t)
    l_grid <- (1:N_l - 0.5) * Delta_l
    a_grid <- (0:N_t) * Delta_t

    # Map observed Length to nearest grid cell index
    length_index <- pmax(1L, pmin(N_l, as.integer(floor(surveys$Length / Delta_l + 0.5))))

    # Unique surveys and indexing
    survey_levels <- sort(unique(surveys$survey_date))
    survey_index <- match(surveys$survey_date, survey_levels)

    # Observations as vectors, aggregated
    obs_df <- data.frame(
        s = survey_index,
        j = length_index,
        K = as.integer(surveys$K),
        count = as.numeric(surveys$count)
    )
    obs_df <- stats::aggregate(count ~ s + j + K, data = obs_df, sum)
    obs_survey_index <- as.integer(obs_df$s)
    obs_length_index <- as.integer(obs_df$j)
    obs_K <- as.integer(obs_df$K)
    obs_count <- as.numeric(obs_df$count)

    # Prepare TMB data and parameters
    tmb_data <- list(
        obs_survey_index = obs_survey_index,
        obs_length_index = obs_length_index,
        obs_K = obs_K,
        obs_count = obs_count,
        survey_dates = as.numeric(survey_levels),
        l_grid = as.numeric(l_grid),
        a_grid = as.numeric(a_grid),
        spawning_mu = as.numeric(pars$spawning_mu),
        spawning_kappa = as.numeric(pars$spawning_kappa),
        annuli_date = as.numeric(pars$annuli_date),
        annuli_min_age = as.numeric(pars$annuli_min_age),
        Delta_l = as.numeric(Delta_l),
        Delta_t = as.numeric(Delta_t),
        log_eps = log(1e-9)
    )

    parameter_names <- c("k", "L_inf", "d", "m")
    tmb_parameters <- pars[parameter_names]

    lower_limit = c(k = 1e-6, L_inf = 1e-3, d = 1e-6, m = 1e-6)
    if (!all(names(lower) %in% parameter_names)) {
        bad_names <- setdiff(names(lower), parameter_names)
        stop("You cannot specify a lower limit on: ", paste(bad_names, collapse = ", "))
    }
    lower_limit[names(lower)] <- lower[names(lower)]

    upper_limit = c(k = Inf, L_inf = Inf, d = Inf, m = Inf)
    if (!all(names(upper) %in% parameter_names)) {
        bad_names <- setdiff(names(upper), parameter_names)
        stop("You cannot specify an upper limit on: ", paste(bad_names, collapse = ", "))
    }
    upper_limit[names(upper)] <- upper[names(upper)]

    obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = tmb_parameters,
        DLL = "growthEstimation",
        silent = TRUE
    )

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

    sdr <- try(TMB::sdreport(obj), silent = TRUE)

    list(
        pars = pars,
        opt = opt,
        sdreport = sdr,
        obj = obj,
        data = tmb_data
    )
}



