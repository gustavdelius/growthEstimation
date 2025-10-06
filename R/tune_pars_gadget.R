#' Interactively tune model parameters with a Shiny gadget
#'
#' Launches a small Shiny gadget with sliders in the left sidebar for each
#' numeric element of `pars` and a live-updating plot of
#' `plot_age_likelihood(pars, age_at_length)` in the main panel.
#'
#' @param pars A named list of scalar numeric parameters. Typical entries are
#'   `k`, `L_inf`, `d`, `m`, `spawning_mu`, `spawning_kappa`, `annuli_date`,
#'   and `annuli_min_age`.
#' @param age_at_length A data frame of observed age-at-length data as used by
#'   `plot_age_likelihood()` (see `Cod_CS_age_at_length`).
#'
#' @return Invisibly returns the final parameters after the gadget is closed.
#'   While the gadget is running it displays the plot interactively.
#'
#' @examples
#' # Not run: requires interactive session
#' # data("Cod_CS_age_at_length")
#' # data("Cod_CS_pars")
#' # tune_pars(Cod_CS_pars, Cod_CS_age_at_length)
#'
#' @export
tune_pars <- function(pars, age_at_length) {

	stopifnot(is.list(pars))

	# Keep only scalar numerics
	numeric_pars <- pars[vapply(pars, function(x) is.numeric(x) && length(x) == 1, logical(1))]
	if (length(numeric_pars) == 0) stop("No scalar numeric parameters found in `pars`.")

	# Helpers for slider rules
	is_fixed_01 <- function(name) grepl("^(annuli_date|annuli_min_age|spawning_mu)$", name)
	is_dynamic  <- function(name) grepl("^(k|L_inf|d|m|r|spawning_kappa)$", name)
	compute_dynamic_range <- function(value_num) {
		minv <- max(0, value_num / 2)
		maxv <- max(1, value_num * 2)
		step <- max((maxv - minv) / 200, 1e-4)
		list(min = minv, max = maxv, step = step)
	}

	# Create a slider per parameter according to the simplified rules
	make_slider <- function(name, value) {
		value_num <- as.numeric(value)
		if (is_fixed_01(name)) {
			minv <- 0; maxv <- 1; step <- 0.001
		} else if (is_dynamic(name)) {
			rng <- compute_dynamic_range(value_num)
			minv <- rng$min; maxv <- rng$max; step <- rng$step
		} else {
			# Default for any other parameter: fixed [0, 1]
			minv <- 0; maxv <- 1; step <- 0.001
		}
		shiny::sliderInput(inputId = name, label = name,
						 min = minv, max = maxv, value = value_num, step = step)
	}

	ui <- shiny::fluidPage(
		shiny::titlePanel("Tune parameters"),
		shiny::sidebarLayout(
			shiny::sidebarPanel(
				shiny::tagList(
					lapply(names(numeric_pars), function(nm) make_slider(nm, numeric_pars[[nm]])),
					shiny::actionButton("done", "Close and return parameters")
				)
			),
			shiny::mainPanel(
				shiny::plotOutput("lik_plot", height = "600px")
			)
		)
	)

	server <- function(input, output, session) {
		current_pars <- shiny::reactive({
			updated <- pars
			for (nm in names(numeric_pars)) {
				updated[[nm]] <- input[[nm]]
			}
			updated
		})

		# Dynamically update ranges for selected parameters when their values change
		for (nm in names(numeric_pars)) {
			if (is_dynamic(nm)) {
				local({
					this_nm <- nm
					shiny::observeEvent(input[[this_nm]], {
						value_num <- as.numeric(input[[this_nm]])
						rng <- compute_dynamic_range(value_num)
						new_value <- min(max(value_num, rng$min), rng$max)
						shiny::updateSliderInput(session, this_nm, min = rng$min, max = rng$max, step = rng$step, value = new_value)
					}, ignoreInit = TRUE)
				})
			}
		}

		output$lik_plot <- shiny::renderPlot({
		    plot_age_likelihood(current_pars(), age_at_length) +
            theme_minimal(base_size = 16)
		})

		# Close gadget and return updated parameters when clicking the button
		shiny::observeEvent(input$done, {
			shiny::stopApp(current_pars())
		})
	}

	# Always run in the system browser
	app <- shiny::shinyApp(ui, server)
	res <- shiny::runApp(app, launch.browser = TRUE)
	invisible(res)
}


