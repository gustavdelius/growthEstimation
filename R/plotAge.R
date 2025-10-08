#' Plot signed log-likelihood contributions from age_at_length with mean K lines
#'
#' Creates a heatmap of the signed contribution of each cell to the NLL,
#' and overlays the mean observed and expected K values.
#' @param pars A list of parameters.
#' @param age_at_length Data frame of age-at-length observations.
#' @return A `ggplot2` object.
#' @export
plot_age_likelihood <- function(pars, age_at_length) {
    # Declare variables to avoid R CMD check warnings
    length <- K <- TotalObserved <- TotalExpected <- SignedNegLogLik <- MeanK <- Source <- NULL

    # Simulate age data
    contributions_df <- get_age_log_likelihood(pars, age_at_length)
    # Convert factors to numeric for plotting
    contributions_df$length <- as.numeric(as.character(contributions_df$length))
    contributions_df$K <- as.numeric(as.character(contributions_df$K))

    # Calculate mean K lines ---
    mean_lines_df <- contributions_df |>
        group_by(length) |>
        summarise(
            # Weighted mean of K by the number of observed/expected fish
            Observed = sum(K * TotalObserved) / sum(TotalObserved),
            Model = sum(K * TotalExpected) / sum(TotalExpected),
            .groups = 'drop'
        ) |>
        # Reshape the data for easy plotting with ggplot
        pivot_longer(
            cols = c("Observed", "Model"),
            names_to = "Source",
            values_to = "MeanK"
        )

    # Determine symmetric color scale limits for the heatmap
    max_abs_val <- max(abs(contributions_df$SignedNegLogLik), na.rm = TRUE)

    # Calculate total negative log likelihood
    total_nll <- sum(contributions_df$TotalNegLogLik, na.rm = TRUE)

    # Create the plot
    p <- ggplot(contributions_df, aes(x = length, y = K)) +
        # Heatmap layer for the misfit
        geom_tile(aes(fill = SignedNegLogLik)) +
        # Line layers for the mean K values
        geom_line(
            data = mean_lines_df,
            aes(y = MeanK, color = Source, linetype = Source),
            linewidth = 1
        ) +
        # --- Define scales and labels ---
        scale_fill_gradient2(
            low = "red",
            mid = "white",
            high = "blue",
            midpoint = 0,
            limit = c(-max_abs_val, max_abs_val),
            name = "Direction & Magnitude\nof Misfit (NLL)"
        ) +
        scale_color_manual(
            name = "Mean K",
            values = c("Observed" = "black", "Model" = "darkgreen")
        ) +
        scale_linetype_manual(
            name = "Mean K",
            values = c("Observed" = "solid", "Model" = "solid")
        ) +
        labs(
            title = "Model Fit Diagnostic",
            subtitle = paste0("Color shows direction (Blue: Obs > Exp, Red: Obs < Exp). Intensity shows magnitude of misfit.\nTotal NLL = ", sprintf("%.2f", total_nll)),
            x = "Fish length (cm)",
            y = "Annuli Count (K)"
        ) +
        theme_minimal()

    return(p)
}

#' Plot observed and predicted length distributions
#'
#' Compares the observed length distribution from data to the model prediction
#' from steady-state density, accounting for size selectivity.
#' @param pars A list containing model parameters to pass to
#'   `solve_pde_steady_state()`: k, L_inf, d, m, r (and optionally vB_min_size),
#'   as well as optional selectivity parameters l50 and l25 (or ratio).
#' @param length_freq A data frame with columns `length` and `count`, containing
#'   observed length frequencies.
#' @param Delta_l Width of size bins (cm). Default is 1.
#' @param l_max The maximum length to consider. If NULL (default), it is set to
#'   110% of the maximum observed length.
#' @return A `ggplot2` object comparing observed and predicted length distributions.
#' @export
#' @examples
#' # Example using the Cod Celtic Sea data
#' data(Cod_CS_pars)
#' data(Cod_CS_length_freq)
#' # Without selectivity
#' plot_length(Cod_CS_pars, Cod_CS_length_freq)
#' # With selectivity
#' pars <- Cod_CS_pars
#' pars$l50 <- 40
#' pars$l25 <- 30
#' plot_length(pars, Cod_CS_length_freq)
plot_length <- function(pars, length_freq, Delta_l = 1, l_max = NULL) {
    # Declare variables to avoid R CMD check warnings
    length <- Source <- Frequency <- NULL

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
    # Calculate l25 from ratio if needed (consistent with TMB implementation)
    if (!is.null(pars$l50)) {
        if (is.null(pars$l25) && !is.null(pars$ratio)) {
            pars$l25 <- pars$ratio * pars$l50
        }

        if (!is.null(pars$l25)) {
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
    }

    # Create data frame for model prediction
    model_df <- data.frame(
        length = l_grid,
        density = u_steady
    )

    # Normalize both observed and predicted to frequencies (proportions)
    total_observed <- sum(length_freq$count)
    total_model <- sum(model_df$density)

    # Create combined data frame for plotting
    observed_df <- data.frame(
        length = length_freq$length,
        Frequency = length_freq$count / total_observed,
        Source = "Observed"
    )

    predicted_df <- data.frame(
        length = model_df$length,
        Frequency = model_df$density / total_model,
        Source = "Model"
    )

    plot_df <- rbind(observed_df, predicted_df)

    # Create the plot
    p <- ggplot(plot_df, aes(x = length, y = Frequency, color = Source, linetype = Source)) +
        geom_line(linewidth = 1) +
        scale_color_manual(
            values = c("Observed" = "black", "Model" = "darkgreen")
        ) +
        scale_linetype_manual(
            values = c("Observed" = "solid", "Model" = "dashed")
        ) +
        labs(
            title = "Length Distribution: Observed vs Model",
            x = "Length (cm)",
            y = "Relative Frequency"
        ) +
        theme_minimal() +
        theme(
            legend.position = "top",
            legend.title = element_blank()
        )

    return(p)
}
