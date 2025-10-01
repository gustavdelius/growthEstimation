#' Fit growth parameters by minimizing the negative log likelihood with TMB
#'
#' Optimizes k, L_inf, d, m, r, l50 and l25/l50 using nlminb().
#' Spawning parameters, annuli_date, and annuli_min_age are treated as data.
#'
#' @param pars List of parameters
#' @param age_at_length Data frame with columns survey_date (numeric), Length, K, count.
#' @param length_freq Data frame with columns survey_date (numeric), Length, and count.
#' @param age_weight Weight given to the age_at_length data
#' @param freq_weight Weight given to the size frequency data
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
    age_at_length,
    length_freq,
    age_weight = 1,
    freq_weight = 1,
    gear_type = 0,
    Delta_l = 1,
    Delta_t = 0.05,
    lower = c(),
    upper = c()
) {
    stopifnot(is.data.frame(age_at_length))
    stopifnot(all(c("survey_date", "Length", "K", "count") %in% names(age_at_length)))
    # If the data was not already aggregated, do it now
    age_at_length |>
        group_by(survey_date, Length, K) |>
        summarise(count = sum(count))

    stopifnot(is.data.frame(length_freq))
    stopifnot(all(c("survey_date", "Length", "count") %in% names(length_freq)))
    # If the data was not already aggregated, do it now
    length_freq |>
        group_by(survey_date, Length) |>
        summarise(count = sum(count))

    # Ensure vB_min_size present in pars
    if (is.null(pars$vB_min_size)) {
        pars$vB_min_size <- as.numeric(min(age_at_length$Length))
    }

    # Build grids similarly to getLogLik()
    l_max_age <- ceiling(max(age_at_length$Length) * 1.1)
    l_max_freq <- ceiling(max(length_freq$Length) * 1.1)
    l_max <- max(l_max_age, l_max_freq)
    t_max <- 100 # set to very high value for the time being because
    # we need the cohort to have very low abundance at final time
    N_l <- ceiling(l_max / Delta_l)
    N_t <- ceiling(t_max / Delta_t)
    # l_grid consists of the centres of the size bins at which we will
    # evaluate the number density
    l_grid <- (1:N_l - 0.5) * Delta_l
    a_grid <- (0:N_t) * Delta_t

    # Map observed Length to nearest grid cell index
    length_index <- pmax(1L, pmin(N_l, as.integer(floor(age_at_length$Length / Delta_l + 0.5))))

    # Unique surveys and indexing for age at length
    survey_levels_age <- sort(unique(age_at_length$survey_date))
    survey_index <- match(age_at_length$survey_date, survey_levels_age)

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

    #Create matrix for length freq data, currently only works with Delta_l = 1
    #cm if Delta_l does not = 1 cm and the survey is not in 1 cm bins
    #then this code will not function properly

    library(dplyr)
    library(tidyr)
    # survey dates and max length
    survey_levels_freq <- sort(unique(length_freq$survey_date))

    # full grid of survey_date x Length
    grid <- expand.grid(
        survey_date = survey_levels_freq,
        Length = 0:(l_max-1)
    )

    # aggregate counts (sum over duplicates)
    survey_counts <- length_freq%>%
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
    rownames(count_matrix) <- survey_levels_freq
    colnames(count_matrix) <- 0:(l_max-1)

    # Prepare TMB data and parameters
    tmb_data <- list(
        obs_survey_index = obs_survey_index,
        obs_length_index = obs_length_index,
        obs_K = obs_K,
        obs_count = obs_count,
        count_matrix = count_matrix,
        age_survey_dates = as.numeric(survey_levels_age),
        freq_survey_dates = as.numeric(survey_levels_freq),
        age_weight = age_weight,
        freq_weight = freq_weight,
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
        gear_type = gear_type
    )

    parameter_names <- c("k", "L_inf", "d", "m", "r", "l50","ratio")
    tmb_parameters <- pars[parameter_names]

    lower_limit = c(k = 1e-6, L_inf = 1e-3, d = 1e-6, m = 1e-6, r = 1e-6,
                    l50 = 1e-6, ratio = 1e-6)
    if (!all(names(lower) %in% parameter_names)) {
        bad_names <- setdiff(names(lower), parameter_names)
        stop("You cannot specify a lower limit on: ", paste(bad_names, collapse = ", "))
    }
    lower_limit[names(lower)] <- lower[names(lower)]

    upper_limit = c(k = Inf, L_inf = Inf, d = Inf, m = Inf, r = Inf,
                    l50 = Inf, ratio = Inf)
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
