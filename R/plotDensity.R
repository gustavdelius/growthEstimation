# Avoid R CMD check notes for non-standard evaluation in ggplot2
utils::globalVariables(c("size", "u", "year"))

#' Plot density function as an interactive 3D surface
#'
#' @param u A numeric matrix containing the density (rows = time, columns =
#'   size).
#' @param Delta_l Size step size used to compute `u`.
#' @param Delta_t Time step size used to compute `u`.
#' @param l_min Size at which the plot should start.
#' @param ... Further arguments passed to `plotly::add_surface()`.
#'
#' @return A Plotly plot
#' @export
#' @examples
#' \dontrun{
#' pars <- list(k = 0.5, L_inf = 100, d = 1, m = 0.2)
#' u <- getDensity(pars, l_max = 100, Delta_l = 1, t_max = 10, Delta_t = 0.1)
#' plot_density_3d(u, Delta_l = 1, Delta_t = 0.1)
#' }
plot_density_3d <- function(u, Delta_l = 1, Delta_t = 0.05, l_min = 5, ...) {
    if (!is.matrix(u)) stop("u must be a numeric matrix.")

    n_time <- nrow(u)
    n_size <- ncol(u)

    # Axes corresponding to cell centres and time steps
    l <- (1:n_size - 0.5) * Delta_l
    t <- (0:(n_time - 1)) * Delta_t

    # Truncate data to include only times >= t_min
    size_indices <- which(l >= l_min)
    if (length(size_indices) == 0) {
        stop("No data available for l >= l_min. Try reducing l_min.")
    }

    # Apply truncation to u and t
    u_truncated <- u[,size_indices , drop = FALSE]
    l_truncated <- l[size_indices]

    # Interactive 3D using plotly
    p <- plotly::plot_ly(x = t, y = l_truncated, z = t(u_truncated))
    p <- plotly::add_surface(p, ...)
    p <- plotly::layout(p, scene = list(xaxis = list(title = "Time t"),
                                        yaxis = list(title = "Size l"),
                                        zaxis = list(title = "u(l, t)")))
    return(p)
}

#' Plot density function snapshots at yearly intervals
#'
#' Uses `ggplot2` to create a 2D line plot of `u(l, t)` for snapshots at yearly
#' intervals, starting at year 1. The matrix returned by `getDensity()` has rows
#' corresponding to time steps (starting at t = 0 in the first row) and columns
#' to size classes.
#'
#' @param u A numeric matrix as returned by `getDensity()` (rows = time,
#'   columns = size).
#' @param Delta_l Size step size used to compute `u`.
#' @param Delta_t Time step size used to compute `u`.
#' @param l_offset Optional offset added to the size axis labels (default 0).
#' @param years Optional numeric vector of years to plot (integers). If not
#'   supplied, uses all integers from 1 up to the largest full year covered by
#'   the simulation time.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' pars <- list(k = 0.5, L_inf = 100, d = 0.01, m = 0.01)
#' u <- getDensity(pars, l_max = 100, Delta_l = 1, t_max = 10, Delta_t = 0.1)
#' plot_density_2d(u, Delta_l = 1, Delta_t = 0.1)
#' }
plot_density_2d <- function(u, Delta_l = 1, Delta_t = 0.05, l_offset = 0, years = NULL) {
    if (!is.matrix(u)) stop("u must be a numeric matrix as returned by getDensity().")

    n_time <- nrow(u)
    n_size <- ncol(u)

    # Axes
    l <- (1:n_size - 0.5) * Delta_l + l_offset
    t <- (0:(n_time - 1)) * Delta_t

    # Determine default years: 1, 2, ..., floor(max(t))
    if (is.null(years)) {
        max_full_year <- floor(max(t))
        years <- seq.int(from = 1L, to = max_full_year, by = 1L)
    }
    if (length(years) == 0) stop("No full years available to plot with the given Delta_t and duration.")

    # For each requested year, find the nearest time index
    nearest_idx <- vapply(years, function(y) {
        which.min(abs(t - y))
    }, integer(1))

    # Build tidy data frame
    u_values <- u[nearest_idx, , drop = FALSE]
    df_list <- lapply(seq_along(years), function(i) {
        data.frame(size = l, u = as.numeric(u_values[i, ]), year = years[i])
    })
    df <- do.call(rbind, df_list)

    # Plot
    p <- ggplot(df, aes(x = size,
                        y = u,
                        colour = factor(year),
                        group = year)) +
        geom_line() +
        labs(x = "Size l", y = "u(l, t)", colour = "Year") +
        theme_minimal()

    return(p)
}
