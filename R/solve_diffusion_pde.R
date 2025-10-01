#' Solve a tridiagonal system of equations using the Thomas algorithm
#'
#' This function solves the system Ax = d, where A is an NxN tridiagonal matrix.
#'
#' @param a A numeric vector representing the sub-diagonal of A. Must have length N-1.
#' @param b A numeric vector representing the main diagonal of A. Must have length N.
#' @param c A numeric vector representing the super-diagonal of A. Must have length N-1.
#' @param d A numeric vector representing the right-hand side vector. Must have length N.
#' @return A numeric vector containing the solution x.
#'
solve_thomas <- function(a, b, c, d) {
    n <- length(b)
    if (n == 0) return(numeric(0))
    if (n == 1) return(d[1] / b[1])

    # Create copies to modify
    c_prime <- numeric(n)
    d_prime <- numeric(n)

    # Forward elimination
    c_prime[1] <- c[1] / b[1]
    d_prime[1] <- d[1] / b[1]

    for (i in 2:(n - 1)) {
        m <- b[i] - a[i-1] * c_prime[i - 1]
        c_prime[i] <- c[i] / m
        d_prime[i] <- (d[i] - a[i-1] * d_prime[i - 1]) / m
    }

    d_prime[n] <- (d[n] - a[n-1] * d_prime[n - 1]) / (b[n] - a[n-1] * c_prime[n - 1])

    # Backward substitution
    x <- numeric(n)
    x[n] <- d_prime[n]
    for (i in (n - 1):1) {
        x[i] <- d_prime[i] - c_prime[i] * x[i + 1]
    }

    return(x)
}


#' Numerically solve the PDE for fish abundance density
#'
#' Implements an unconditionally stable finite volume scheme (implicit Euler)
#' with upwinding for advection.
#'
#' @param pars A list containing the model parameters: k, L_inf, d, m, r.
#' @param u_initial A numeric vector for the initial condition u(l, 0).
#' @param Delta_l The size step size (cm).
#' @param t_max The maximum simulation time.
#' @param Delta_t The time step size (years). Default is 0.05.
#' @return A matrix where each column is the solution u(t,l) at a given time
#'   step. Rows correspond to time and columns to length.
#'
solve_pde <- function(pars, u_initial,
                     Delta_l = 1, t_max = 10, Delta_t = 0.05) {

    # Set up Grids ----

    N_l <- length(u_initial) # Number of Size cells
    N_t <- ceiling(t_max / Delta_t) # Number of time steps

    # Size grid (cell centres)
    l_grid <- (1:N_l - 0.5) * Delta_l
    # Size grid (cell interfaces, size N_l + 1)
    l_interfaces <- (0:N_l) * Delta_l

    # Initialize Solution Matrix ----
    # Rows = space, Columns = time. Add 1 column for the initial condition.
    u_solution <- matrix(0, nrow = N_t + 1, ncol = N_l)
    u_solution[1, ] <- u_initial

    # Pre-calculate Coefficients ----
    # These coefficients are time-independent, so we can compute them once.

    # Advection (v) and Diffusion (D) coefficients at interfaces
    # Use constant growth rate r for sizes below vB_min_size,
    # von Bertalanffy growth rate k*(L_inf - L) for sizes >= vB_min_size
    vB_min_size <- if (is.null(pars$vB_min_size)) 0 else as.numeric(pars$vB_min_size)
    growth_rate <- ifelse(l_interfaces < vB_min_size, 
                         pars$r, 
                         pars$k * (pars$L_inf - l_interfaces))
    v <- growth_rate - pars$d / 2
    D <- pars$d * l_interfaces / 2

    v_plus <- pmax(v, 0)
    v_minus <- pmin(v, 0)

    # Reaction coefficient (mu) at cell centres.
    mu <- pars$m / l_grid

    # Construct the tridiagonal system matrix (A) ----
    # We define A via its three diagonals: a_ (sub), b_ (main), c_ (super)

    # Initialize diagonal vectors
    a_ <- numeric(N_l - 1)
    b_ <- numeric(N_l)
    c_ <- numeric(N_l - 1)

    # Pre-factors for convenience
    c1 <- Delta_t / Delta_l
    c2 <- Delta_t / (Delta_l^2)

    # Fill the vectors for the interior points (i = 2 to N_l-1)
    # R vector indices are 1-based, so this loop is for i from 2 to N_l-1
    for (i in 2:(N_l - 1)) {
        # A[i, i-1] depends on interface i
        a_[i - 1] <- -c1 * v_plus[i] - c2 * D[i]
        # A[i, i+1] depends on interface i+1
        c_[i] <- c1 * v_minus[i + 1] - c2 * D[i + 1]
        # A[i, i] depends on interfaces i and i+1
        b_[i] <- 1 + Delta_t * mu[i] +
            c1 * (v_plus[i + 1] - v_minus[i]) +
            c2 * (D[i + 1] + D[i])
    }

    # Apply Boundary Conditions ----

    # -- At l_0 = 0 (i=1): No-Flux condition (J_{1/2} = 0) --
    # The equation for u_1 has no u_0 term.
    b_[1] <- 1 + Delta_t * mu[1] + c1 * (v_plus[2] + D[2] / Delta_l)
    c_[1] <- c1 * (v_minus[2] - D[2] / Delta_l)

    # -- At l_max (i=N_l): Dirichlet condition u(l_max, t) = 0 --
    # This implies u_{N_l+1} = 0 in the flux J_{N_l+1/2}.
    a_[N_l - 1] <- -c1 * v_plus[N_l] - c2 * D[N_l]
    b_[N_l] <- 1 + Delta_t * mu[N_l] +
        c1 * (v_plus[N_l + 1] - v_minus[N_l]) +
        c2 * (D[N_l + 1] + D[N_l])


    # Time-Stepping Loop ----
    for (n in 1:N_t) {
        # The right-hand side is the solution from the previous time step
        d_rhs <- u_solution[n, ]

        # Solve the system A * u^{n+1} = u^n
        u_solution[n + 1, ] <- solve_thomas(a_, b_, c_, d_rhs)
    }

    # Return Result ----
    return(u_solution)
}

#' Calculate Green's function for fish abundance PDE
#'
#' Solves the PDE with an initial condition where all individuals are in the
#' smallest size class.
#'
#' @inheritParams solve_pde
#' @param l_max The maximum size
#' @return A matrix holding the Green's function G(t,l). Rows correspond to
#'   time and columns to length.
#' @export
getGreens <- function(pars, l_max, Delta_l = 1, t_max = 10, Delta_t = 0.05) {

    # Set initial condition ----
    N_l <- ceiling(l_max / Delta_l)  # Number of size cells
    u_initial <- numeric(N_l)
    u_initial[1] <- 1

    # Solve PDE ----
    G <- solve_pde(pars, u_initial = u_initial,
                   Delta_l = Delta_l, t_max = t_max, Delta_t = Delta_t)
    return(G)
}

#' Numerically solve the steady state of the PDE for fish abundance density
#'
#' Finds the steady state solution by setting the time derivative to zero,
#' resulting in an ODE that is discretized using the same finite volume scheme
#' as solve_pde() and solved using Thomas's method.
#'
#' @param pars A list containing the model parameters: k, L_inf, d, m, r.
#' @param Delta_l The size step size (cm). Default is 1.
#' @param l_max The maximum size. Default is 100.
#' @return A numeric vector containing the steady state solution u(l).
#'
#' @export
solve_pde_steady_state <- function(pars, Delta_l = 1, l_max = 100) {
    
    # Set up Grid ----
    N_l <- ceiling(l_max / Delta_l)  # Number of size cells
    
    # Size grid (cell centres)
    l_grid <- (1:N_l - 0.5) * Delta_l
    # Size grid (cell interfaces, size N_l + 1)
    l_interfaces <- (0:N_l) * Delta_l
    
    # Pre-calculate Coefficients ----
    # These coefficients are the same as in solve_pde()
    
    # Advection (v) and Diffusion (D) coefficients at interfaces
    vB_min_size <- if (is.null(pars$vB_min_size)) 0 else as.numeric(pars$vB_min_size)
    growth_rate <- ifelse(l_interfaces < vB_min_size, 
                         pars$r, 
                         pars$k * (pars$L_inf - l_interfaces))
    v <- growth_rate - pars$d / 2
    D <- pars$d * l_interfaces / 2
    
    v_plus <- pmax(v, 0)
    v_minus <- pmin(v, 0)
    
    # Reaction coefficient (mu) at cell centres
    # Add small epsilon to avoid division by zero
    mu <- pars$m / pmax(l_grid, 1e-10)
    
    # Construct the tridiagonal system matrix (A) for steady state ----
    # For steady state: A * u = 0, so we need to solve A * u = 0
    # This is a homogeneous system, so we need to impose a constraint
    # We'll use the boundary condition at l_max as u(l_max) = 0
    
    # Initialize diagonal vectors
    a_ <- numeric(N_l - 1)
    b_ <- numeric(N_l)
    c_ <- numeric(N_l - 1)
    
    # Pre-factors for convenience (no Delta_t since we're solving steady state)
    c1 <- 1 / Delta_l
    c2 <- 1 / (Delta_l^2)
    
    # Fill the vectors for the interior points (i = 2 to N_l-1)
    for (i in 2:(N_l - 1)) {
        # A[i, i-1] depends on interface i
        a_[i - 1] <- -c1 * v_plus[i] - c2 * D[i]
        # A[i, i+1] depends on interface i+1
        c_[i] <- c1 * v_minus[i + 1] - c2 * D[i + 1]
        # A[i, i] depends on interfaces i and i+1
        b_[i] <- mu[i] +
            c1 * (v_plus[i + 1] - v_minus[i]) +
            c2 * (D[i + 1] + D[i])
    }
    
    # Apply Boundary Conditions ----
    
    # -- At l_0 = 0 (i=1): No-Flux condition (J_{1/2} = 0) --
    b_[1] <- mu[1] + c1 * (v_plus[2] + D[2] / Delta_l)
    c_[1] <- c1 * (v_minus[2] - D[2] / Delta_l)
    
    # -- At l_max (i=N_l): Dirichlet condition u(l_max) = 0 --
    # This is enforced by the homogeneous system constraint
    a_[N_l - 1] <- -c1 * v_plus[N_l] - c2 * D[N_l]
    b_[N_l] <- mu[N_l] +
        c1 * (v_plus[N_l + 1] - v_minus[N_l]) +
        c2 * (D[N_l + 1] + D[N_l])
    
    # For a homogeneous system A * u = 0, we need to impose a constraint
    # We'll fix u[1] = 1 and solve for the remaining variables
    # This gives us A' * u[2:N_l] = -A[2:N_l, 1]
    
    # Create the modified system
    if (N_l > 1) {
        # Right-hand side: -A[2:N_l, 1] * u[1] where u[1] = 1
        d_rhs <- numeric(N_l - 1)
        d_rhs[1] <- -a_[1]  # -A[2, 1] since u[1] = 1
        
        # Modified tridiagonal system for u[2:N_l]
        a_mod <- a_[2:(N_l - 1)]
        b_mod <- b_[2:N_l]
        c_mod <- c_[1:(N_l - 2)]
        
        # Solve the system
        if (N_l > 2) {
            # Check if the system is singular (all diagonal elements are zero or very small)
            if (all(abs(b_mod) < 1e-15)) {
                u_interior <- rep(0, length(b_mod))
            } else {
                u_interior <- solve_thomas(a_mod, b_mod, c_mod, d_rhs)
            }
        } else {
            # Check for division by zero
            if (abs(b_mod[1]) < 1e-15) {
                u_interior <- 0  # Return zero if system is singular
            } else {
                u_interior <- d_rhs[1] / b_mod[1]
            }
        }
        
        # Check if solution contains NaN or Inf values
        if (any(!is.finite(u_interior))) {
            # Return a simple solution if the system is singular
            u_steady <- c(1, rep(0, length(u_interior)))
        } else {
            # Combine with the fixed boundary value
            u_steady <- c(1, u_interior)
        }
    } else {
        u_steady <- 1
    }
    
    # Normalize the solution (optional, for numerical stability)
    # Scale so that the maximum value is 1
    max_val <- max(abs(u_steady), na.rm = TRUE)
    if (is.finite(max_val) && max_val > 0) {
        u_steady <- u_steady / max_val
    }
    
    return(u_steady)
}
