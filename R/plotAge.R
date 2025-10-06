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
