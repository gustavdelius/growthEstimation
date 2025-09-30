#' Fit growth parameters by minimizing the negative log likelihood with TMB
#'
#' Optimizes k, L_inf, d, m, l50 and ratio using nlminb().
#' Spawning parameters are treated as data.
#' @param pars List of parameters with k, L_inf, d, m, spawning_mu, and spawning_kappa
#' @param surveys Data frame with columns survey_date (numeric), Length, and count.
#' @param gear_type Sigmoidal (0) and knife edge (1) are the only supported types, default 0.
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
    stopifnot(all(c("survey_date", "Length", "count") %in% names(surveys)))
    surveys <- as.data.frame(surveys)

    # Build grids similarly to getLogLik()
    l_max <- ceiling(max(surveys$Length) * 1.1)
    t_max <- 100
    N_l <- ceiling(l_max / Delta_l)
    N_t <- ceiling(t_max / Delta_t)
    l_grid <- (1:N_l - 0.5) * Delta_l
    a_grid <- (0:N_t) * Delta_t

    #matix
    library(dplyr)
    library(tidyr)

    # survey dates and max length
    survey_levels <- sort(unique(surveys$survey_date))

    # full grid of survey_date x Length
    grid <- expand.grid(
        survey_date = survey_levels,
        Length = 0:(l_max-1)
    )

    # aggregate counts (sum over duplicates)
    survey_counts <- surveys %>%
        group_by(survey_date, Length) %>%
        summarise(count = sum(count), .groups = "drop")

    # merge counts into full grid, fill missing with 0
    survey_grid <- grid %>%
        left_join(survey_counts, by = c("survey_date", "Length")) %>%
        mutate(count = replace_na(count, 0))

    # reshape to wide matrix
    count_matrix <- survey_grid %>%
        pivot_wider(
            names_from = Length,
            values_from = count,
            values_fill = 0
        )

    # convert to plain matrix with proper dimnames
    count_matrix <- as.matrix(count_matrix[-1])
    rownames(count_matrix) <- survey_levels
    colnames(count_matrix) <- 0:(l_max-1)

    # Prepare TMB data and parameters
    tmb_data <- list(
        count_matrix = count_matrix,
        survey_dates = as.numeric(survey_levels),
        l_grid = as.numeric(l_grid),
        a_grid = as.numeric(a_grid),
        spawning_mu = as.numeric(pars$spawning_mu),
        spawning_kappa = as.numeric(pars$spawning_kappa),
        annuli_date = as.numeric(pars$annuli_date),
        Delta_l = as.numeric(Delta_l),
        Delta_t = as.numeric(Delta_t),
        log_eps = log(1e-9),
        gear_type = gear_type
    )

    parameter_names <- c("k", "L_inf", "d", "m", "l50", "ratio")
    tmb_parameters <- pars[parameter_names]

    lower_limit = c(k = 1e-6, L_inf = 1e-3, d = 1e-6, m = 1e-6, l50 = 1e-3, ratio = 1e-3)
    if (!all(names(lower) %in% parameter_names)) {
        bad_names <- setdiff(names(lower), parameter_names)
        stop("You cannot specify a lower limit on: ", paste(bad_names, collapse = ", "))
    }
    lower_limit[names(lower)] <- lower[names(lower)]

    upper_limit = c(k = Inf, L_inf = Inf, d = Inf, m = Inf, l50 = Inf, ratio = Inf)
    if (!all(names(upper) %in% parameter_names)) {
        bad_names <- setdiff(names(upper), parameter_names)
        stop("You cannot specify an upper limit on: ", paste(bad_names, collapse = ", "))
    }
    upper_limit[names(upper)] <- upper[names(upper)]

    obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = tmb_parameters,
        DLL = "length_freq",
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
