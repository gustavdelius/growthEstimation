#' Get log likelihood of age observations
#'
#' @param pars A list containing the model parameters: k, L_inf, d, m, r.
#' @param surveys A data frame with survey age-at-length observations with
#'   columns `survey_date`, `Length`, `K`, and `count`.
#' @param Delta_l Width of size bins (cm). Default is 1.
#' @param Delta_t Time step for the model simulation (years). Default is 0.05.
#'
#' @return A data frame with, for each observed Length-K bin in each survey, the
#'   observed count, expected count under the model, model probability, sample
#'   size, negative log-likelihood contribution, and signed negative
#'   log-likelihood contribution.
#' @export
getLogLik <- function(pars, surveys, Delta_l = 1, Delta_t = 0.05) {
    # Set up grids ----
    l_max <- ceiling(max(surveys$Length) * 1.1) # A bit larger than maximum observed size
    t_max <- max(surveys$K) + 2 # TODO: set t_max only as high as needed
    N_l <- ceiling(l_max / Delta_l) # Number of time steps
    N_t <- ceiling(t_max / Delta_t) # Number of time steps
    l_grid <- (1:N_l) * Delta_l # Size grid (bin start values)
    a_grid <- (0:N_t) * Delta_t   # Age grid (including age = 0)

    # Ensure vB_min_size in pars once for downstream calls
    if (is.null(pars$vB_min_size)) {
        pars$vB_min_size <- as.numeric(min(surveys$Length))
    }
    G <- getGreens(pars, l_max = l_max, Delta_l = Delta_l,
                   t_max = t_max, Delta_t = Delta_t)

    calculate_and_aggregate_likelihood(
        surveys, G = G, a = a_grid, l = l_grid,
        mu = pars$spawning_mu, kappa = pars$spawning_kappa,
        annuli_date = pars$annuli_date, annuli_min_age = pars$annuli_min_age)
}
