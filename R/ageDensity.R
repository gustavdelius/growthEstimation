#' Age-to-Annuli Mapping Function
#'
#' Deterministically maps true age to the expected annuli count, given a survey
#' date, an annual ring-formation day, and a minimum age threshold.
#' @param age_in_years A numeric vector of true ages in years.
#' @param survey_date The survey date as a numeric year (e.g., 2023.25).
#' @param annuli_date Ring formation day as fraction of a year in \[0, 1).
#' @param annuli_min_age The minimum age (years) required to form the first
#'   ring.
#' @return An integer vector with the calculated number of rings (K) for each
#'   age.
#' @export
age_to_annuli <- function(age_in_years, survey_date,
                          annuli_date, annuli_min_age) {
    sapply(age_in_years, function(age) {
        if (age < 0) return(0)
        birth_date <- survey_date - age
        birth_year <- floor(birth_date)
        next_ring_date <- birth_year + annuli_date
        if (next_ring_date <= birth_date) {
            next_ring_date <- (birth_year + 1) + annuli_date
        }
        k_count <- 0
        while (next_ring_date < survey_date) {
            age_at_ring_date <- next_ring_date - birth_date
            if (age_at_ring_date >= annuli_min_age) {
                k_count <- k_count + 1
            }
            next_ring_date <- next_ring_date + 1
        }
        return(k_count)
    })
}

#' Generate model predictions for annuli at length at a specific survey date
#'
#' Calculates the probability P(K | length) that a fish in a given length class
#' is observed with K annuli.
#' @param survey_date Numeric survey date (e.g., 2023.25).
#' @param G Greens function matrix (rows = ages, cols = length classes).
#' @param a Numeric vector of high-resolution ages corresponding to rows of `G`.
#' @param l Numeric vector of lengths corresponding to columns of `G`.
#' @param mu Mean spawning date as a fraction of a year in \[0, 1).
#' @param kappa Spawning concentration parameter.
#' @param annuli_date Ring formation day as fraction of a year in \[0, 1).
#' @param annuli_min_age Minimum age (years) at which the first ring can form.
#' @return A matrix of probabilities P(K|l) with rows named by `length` and
#'   columns by `K`.
#' @export
#' @examples
#' # Minimal schematic example (using toy inputs)
#' a <- seq(0, 3, length.out = 5)
#' l <- seq(10, 30, length.out = 3)
#' G <- matrix(abs(sin(outer(a, l, "+"))), nrow = length(a))
#' annuli_predictions_for_date(2023.5, G, a, l, mu = 0.5, kappa = 3,
#'                                     annuli_date = 0.25, annuli_min_age = 0.5)
annuli_predictions_for_date <- function(
        survey_date, G, a, l, mu, kappa, annuli_date, annuli_min_age) {
    # Population Convolution
    birth_dates <- survey_date - a
    spawning_weights <- get_spawning_density(birth_dates, mu, kappa)
    N_pop <- diag(spawning_weights) %*% G

    # Observation Convolution
    k_for_each_age <- age_to_annuli(a, survey_date, annuli_date, annuli_min_age)
    max_K <- max(k_for_each_age)
    k_bins <- 0:max_K
    N_model <- matrix(0, nrow = length(l), ncol = length(k_bins))
    dimnames(N_model) <- list(length = l, K = k_bins)

    for (k_val in k_bins) {
        age_indices <- which(k_for_each_age == k_val)
        if (length(age_indices) > 0) {
            pop_subset <- N_pop[age_indices, , drop = FALSE]
            N_model[, k_val + 1] <- colSums(pop_subset)
        }
    }

    # Step D: Convert to proportions
    P_model_K_given_l <- prop.table(N_model, margin = 1)
    P_model_K_given_l[is.nan(P_model_K_given_l)] <- 0
    return(P_model_K_given_l)
}

#' Calculate and Aggregate Signed Log-Likelihood Contributions
#'
#' Loops through age_at_length, calculates signed NLL for each, and aggregates
#' the results.
#' @param age_at_length A data frame with survey age-at-length observations with
#'   columns `survey_date`, `length`, `K`, and `count`.
#' @inheritParams annuli_predictions_for_date
#' @return A data frame containing, for each observed length-K bin in each
#'   survey, the observed count, expected count under the model, model
#'   probability, sample size, negative log-likelihood contribution, and signed
#'   negative log-likelihood contribution.
#' @export
age_likelihood <- function(age_at_length, G, a, l, mu, kappa,
                           annuli_date, annuli_min_age) {
    # Declare variables to avoid R CMD check warnings
    length <- K <- N <- Prob <- Expected <- NegLogLik <- SignedNegLogLik <-
        TotalObserved <- TotalExpected <- TotalNegLogLik <- count <- NULL

    # Split the data frame by unique survey date
    age_at_length <- split(age_at_length, age_at_length$survey_date)

    log_lik_contributions <- list()

    for (survey_date_str in names(age_at_length)) {
        survey_date_current <- as.numeric(survey_date_str)

        current_obs_df <- age_at_length[[survey_date_str]]

        # 1. Generate model predictions for this date
        P_model <- annuli_predictions_for_date(
            survey_date_current, G, a, l, mu, kappa,
            annuli_date, annuli_min_age
        )

        # 2. Get sample sizes per length for this survey

        sample_sizes <- current_obs_df |>
            group_by(length) |>
            summarise(N = sum(count, na.rm = TRUE), .groups = 'drop')

        # 3. Convert model probabilities to a long data frame for joining
        P_model_df <- melt(P_model, value.name = "Prob")

        # 4. Join all data together
        likelihood_df <- current_obs_df |>
            left_join(sample_sizes, by = "length") |>
            left_join(P_model_df, by = c("length", "K"))

        # 5. Calculate the signed negative log-likelihood contribution
        epsilon <- 1e-9 # To prevent log(0)
        likelihood_df <- likelihood_df |>
            mutate(
                Expected = N * Prob,
                NegLogLik = -(count * log(Prob + epsilon)),
                SignedNegLogLik = sign(count - Expected) * NegLogLik
            )

        log_lik_contributions[[survey_date_str]] <- likelihood_df
    }

    # Aggregate contributions across all age_at_length surveys
    all_contributions_df <- do.call(rbind, log_lik_contributions)

    total_contributions <- all_contributions_df |>
        group_by(length, K) |>
        summarise(
            TotalObserved = sum(count, na.rm = TRUE),
            TotalExpected = sum(Expected, na.rm = TRUE),
            TotalNegLogLik = sum(NegLogLik, na.rm = TRUE),
            .groups = 'drop'
        ) |>
        mutate(
            # The final signed NLL is the total misfit, with the sign determined
            # by the overall difference between observed and expected counts.
            SignedNegLogLik = sign(TotalObserved - TotalExpected) *
                TotalNegLogLik
        )

    return(total_contributions)
}

#' Simulate a sample from model predictions
#'
#' Draws K values using multinomial sampling with probabilities P(K | length)
#' for each length observed in a given survey.
#' @param P_model_K_given_l Matrix of predicted probabilities P(K | length).
#' @param survey_obs A data frame for one survey with at least a `length`
#'   column.
#' @return An integer vector of simulated K values aligned with `survey_obs`
#'   rows.
#' @export
#' @examples
#' set.seed(1)
#' P <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
#' rownames(P) <- c("20", "30"); colnames(P) <- c("0", "1")
#' survey_obs <- data.frame(length = c(20, 20, 30, 30, 30))
#' simulate_age_sample_from_model(P, survey_obs)
simulate_age_sample_from_model <- function(P_model_K_given_l, survey_obs) {
    simulated_K <- integer(nrow(survey_obs))

    # Get the unique length classes that were actually sampled in this survey
    unique_lengths <- unique(survey_obs$length)

    # Get the possible K values from the column names of the proportion matrix.
    k_values <- as.numeric(colnames(P_model_K_given_l))

    # For each unique length class...
    for (len_val in unique_lengths) {
        indices_to_fill <- which(survey_obs$length == len_val)
        n_fish <- length(indices_to_fill)
        len_char <- as.character(len_val)

        # Get the model's predicted proportions of K for this length
        k_proportions <- P_model_K_given_l[len_char, , drop = TRUE]

        # Check if the length class exists in the model predictions and has valid probabilities
        if (len_char %in% rownames(P_model_K_given_l) && sum(k_proportions, na.rm = TRUE) > 0) {

            # Perform a multinomial random sample
            sim_counts <- rmultinom(1, size = n_fish, prob = k_proportions)

            # Create a vector of the simulated K values
            simulated_k_vector <- rep(k_values, sim_counts)

            # Assign these simulated K's to the correct rows in the output.
            # The length of simulated_k_vector is guaranteed to match n_fish.
            if (length(simulated_k_vector) > 0) {
                simulated_K[indices_to_fill] <- simulated_k_vector
            }
        }
    }
    return(simulated_K)
}
