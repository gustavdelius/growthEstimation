# Here we simulate age-at-size data and see how well we can
# estimate the parameters from the simulated data.

true_pars <- list(
    k = 0.2, L_inf = 90, d = 0.5, m = 4, r = 5,
    l50 = 20, l25 = 15,
    spawning_mu = 0.5, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
)
start_pars <- list(
    k = 0.1, L_inf = 150, d = 0.1, m = 10, r = 1,
    l50 = 10, l25 = 9,
    spawning_mu = 0.5, spawning_kappa = 3,
    annuli_date = 0.25, annuli_min_age = 0.5
)

survey_dates <- c(2023.25, 2023.5, 2023.75)
lengths <- seq(10, 100, by = 1)
n_per_length = 20
n_total <- 5000

set.seed(40)
res <- fit_on_simulated_data(
    true_pars = true_pars,
    start_pars = start_pars,
    survey_dates = survey_dates,
    lengths = lengths,
    n_per_length = n_per_length,
    n_length_total = n_total,
    weight_length = 6,
    lower = c()
)
res$fit$opt$convergence
res$comparison

# ---- Diagnose parameter correlation issues ----

# 2. Extract Hessian and compute correlation matrix
cat("=== Parameter Correlation Analysis ===\n")
hessian <- res$fit$obj$he()
param_names <- names(res$fit$obj$par)

# Check if Hessian is positive definite
hess_eigen <- eigen(hessian, symmetric = TRUE, only.values = TRUE)
cat("Hessian eigenvalues (should all be positive):\n")
print(hess_eigen$values)
cat("\nCondition number (ratio of max to min eigenvalue):",
    max(hess_eigen$values) / min(hess_eigen$values), "\n")
cat("(High condition number >1000 suggests ill-conditioning)\n\n")

# Compute covariance and correlation matrices
if (all(hess_eigen$values > 0)) {
    cov_matrix <- solve(hessian)
    std_errors <- sqrt(diag(cov_matrix))
    cor_matrix <- cov_matrix / outer(std_errors, std_errors)
    rownames(cor_matrix) <- colnames(cor_matrix) <- param_names

    cat("Parameter estimates with standard errors:\n")
    print(data.frame(
        Estimate = res$fit$opt$par,
        Std_Error = std_errors,
        CV = std_errors / abs(res$fit$opt$par)
    ))
    cat("\n")

    cat("Correlation matrix:\n")
    print(round(cor_matrix, 3))
    cat("\n")

    # Flag high correlations
    high_cor_threshold <- 0.95
    high_cors <- which(abs(cor_matrix) > high_cor_threshold &
                       abs(cor_matrix) < 1, arr.ind = TRUE)
    if (length(high_cors) > 0 && nrow(high_cors) > 0) {
        # Remove duplicates (upper triangle only)
        upper_tri_mask <- high_cors[,1] < high_cors[,2]
        high_cors <- high_cors[upper_tri_mask, , drop = FALSE]

        if (nrow(high_cors) > 0) {
            cat("WARNING: High correlations (|r| >", high_cor_threshold, ") detected:\n")
            for (i in 1:nrow(high_cors)) {
                par1 <- param_names[high_cors[i,1]]
                par2 <- param_names[high_cors[i,2]]
                cor_val <- cor_matrix[high_cors[i,1], high_cors[i,2]]
                cat(sprintf("  %s <-> %s: %.3f\n", par1, par2, cor_val))
            }
        } else {
            cat("No extreme correlations detected.\n")
        }
    } else {
        cat("No extreme correlations detected.\n")
    }

    # Visualize correlation matrix
    if (requireNamespace("corrplot", quietly = TRUE)) {
        cat("\nGenerating correlation plot...\n")
        corrplot::corrplot(cor_matrix, method = "color", type = "upper",
                          addCoef.col = "black", tl.col = "black",
                          tl.srt = 45, number.cex = 0.7,
                          title = "Parameter Correlation Matrix",
                          mar = c(0,0,2,0))
    }

} else {
    cat("WARNING: Hessian is not positive definite!\n")
    cat("This indicates the optimization may not have reached a true minimum.\n")
    cat("Negative eigenvalues:", hess_eigen$values[hess_eigen$values <= 0], "\n\n")
}

# 3. Check gradient at optimum
cat("=== Gradient Check ===\n")
gradient <- res$fit$obj$gr(res$fit$opt$par)
cat("Gradient at optimum (should be close to 0):\n")
print(setNames(gradient, param_names))
cat("Max absolute gradient:", max(abs(gradient)), "\n\n")

cat("===================================\n\n")

plot_age_likelihood(res$fit$pars, res$simulated_age_data)
plot_age_likelihood(true_pars, res$simulated_age_data)
plot_length(res$fit$pars, res$simulated_length_data)
plot_length(true_pars, res$simulated_length_data)

u <- get_number_density(res$fit$pars, l_max = 100)
plot_density_2d(u)
plot_density_3d(u, l_min = 10)


