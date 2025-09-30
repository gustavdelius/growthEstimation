#' Example age_at_length data for Cod in the Celtic Sea.
#'
#' @format A data frame with 329 rows and 4 columns. The columns ae:
#' \describe{
#'   \item{survey_date}{The within-year date of the survey (decimal year)}
#'   \item{Length}{Length (cm)}
#'   \item{K}{Number of observed annuli}
#'   \item{count}{Number of fish observed with given length and annuli}
#' }
#' @name Cod_CS_age_at_length
#' @docType data
#' @keywords datasets
"Cod_CS_age_at_length"

#' Example parameter set for Cod in the Celtic Sea.
#'
#' @format A named list of scalar numeric parameters used by the growth and
#' age-at-length likelihood functions. Typical entries are:
#' \describe{
#'   \item{k}{Growth rate parameter for von Bertalanffy growth}
#'   \item{L_inf}{Asymptotic length}
#'   \item{d}{Diffusion/variability parameter}
#'   \item{m}{Mortality parameter}
#'   \item{r}{Constant growth rate parameter for small sizes}
#'   \item{vB_min_size}{Size at which growth becomes von Bertalanffy (cm). Used to
#'     switch between constant growth (below) and von Bertalanffy growth (above).}
#'   \item{spawning_mu}{Mean spawning date within the year (decimal year [0,1))}
#'   \item{spawning_kappa}{Concentration of spawning timing (Von Mises kappa)}
#'   \item{annuli_date}{Within-year date when annuli are laid (decimal year)}
#'   \item{annuli_min_age}{Minimum age (in years) before first annulus is counted}
#' }
#' @name Cod_CS_pars
#' @docType data
#' @keywords datasets
"Cod_CS_pars"
