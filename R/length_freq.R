#' Fit growth parameters by minimizing the negative log likelihood with TMB
#'
#' Optimizes k, L_inf, d, m, and selectivity using nlminb().
#' Spawning parameters are treated as data.
#'
#' @param pars List of parameters with k, L_inf, d, m, spawning_mu, and spawning_kappa
#' @param gear_type Sigmoidal (0) and knife edge (1) are the only supported types, default 0.
#' @param surveys Data frame with columns survey_date (numeric), Length, and count.
#' @param Delta_l Numeric size step (cm), default 1.
#' @param Delta_t Numeric time step (years), default 0.05.
#' @param lower Named numeric vector of lower bounds.
#' @param upper Named numeric vector of upper bounds.
#' @return A list with updated `pars`, optimizer result, sdreport, obj, and data
#'   used.
#' @export
fit_tmb_length_freq <- function(
        pars,
        surveys,
        gear_type= 0,
        Delta_l = 1,
        Delta_t = 0.05,
        lower = c(),
        upper = c()
) {length_freq <- as.data.frame(surveys)

# Build grids similarly to getLogLik()
l_max <- ceiling(max(surveys$Length) * 1.1)
t_max <- 100 ##set to very high value for the time being
N_l <- ceiling(l_max / Delta_l)
N_t <- ceiling(t_max / Delta_t)
l_grid <- (1:N_l - 0.5) * Delta_l
a_grid <- (0:N_t) * Delta_t

# Map observed Length to nearest grid cell index
length_index <- pmax(1L, pmin(N_l, as.integer(floor(surveys$Length / Delta_l + 0.5))))

# Unique surveys and indexing
survey_levels <- sort(unique(surveys$survey_date))
survey_index <- match(surveys$survey_date, survey_levels)

# Function to build zero-padded count_matrix
build_count_matrix <- function(surveys, l_grid) {
    survey_levels <- sort(unique(surveys$survey_date))
    nSurvey <- length(survey_levels)

    # Map survey dates to indices (0-based for TMB)
    obs_survey_index <- match(surveys$survey_date, survey_levels) - 1L

    # Map lengths to shared grid
    obs_length_index <- as.integer(floor(surveys$Length / Delta_l + 0.5) - 1L)

    # Build zero-padded count matrix
    count_matrix <- matrix(0, nrow = nSurvey, ncol = length(l_grid))
    for (i in seq_len(nrow(surveys))) {
        s <- obs_survey_index[i] + 1L
        l <- obs_length_index[i] + 1L
        count_matrix[s, l] <- count_matrix[s, l] + surveys$count[i]
    }
    storage.mode(count_matrix) <- "double"

    # Prepare observation vectors for TMB
    obs_df <- data.frame(
        s = obs_survey_index + 1L,
        l = obs_length_index + 1L,
        count = surveys$count
    )
    obs_df <- stats::aggregate(count ~ s + l, data = obs_df, sum)

    list(
        count_matrix = count_matrix,
        obs_survey_index = as.integer(obs_df$s - 1L),
        obs_length_index = as.integer(obs_df$l - 1L),
        obs_count = as.numeric(obs_df$count),
        survey_levels = survey_levels
    )
}

lf_data <- build_count_matrix(length_freq, l_grid)

# Prepare TMB data and parameters
tmb_data <- list(
    obs_survey_index = lf_data$obs_survey_index,
    obs_length_index = lf_data$obs_length_index,
    obs_count = lf_data$obs_count,
    count_matrix = lf_data$count_matrix,
    survey_dates = as.numeric(lf_data$survey_levels),
    l_grid = as.numeric(l_grid),
    a_grid = as.numeric(a_grid),
    spawning_mu = as.numeric(pars["spawning_mu"]),
    spawning_kappa = as.numeric(pars["spawning_kappa"]),
    Delta_l = as.numeric(Delta_l),
    Delta_t = as.numeric(Delta_t),
    log_eps = log(1e-9),
    gear_type = gear_type
)

parameter_names <- c("k", "L_inf", "d", "m", "l50","ratio")
tmb_parameters <- pars[parameter_names]

lower_limit = c(k = 1e-6, L_inf = 1e-3, d = 1e-6, m = 1, l50 = 1e-6, ratio = 1e-6)
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
    DLL = "Length_freq",
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
