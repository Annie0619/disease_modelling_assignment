# R/fit_nb.R
# -----------------------------------------------------------------------------
# Negative binomial fit of the model to the data (parallel multi-start).
# This file exposes a single user-facing function:
#   fit_nb_to_cases(...)
#
# Nothing here mutates global state. All inputs are passed explicitly.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(future)
  library(furrr)
})

# -----------------------------------------------------------------------------
# Public entry point (placeholder for now; we will fill in later steps)
# -----------------------------------------------------------------------------
# Arguments:
#   mdl                : model object from build_model()
#   y0_burn            : state vector at start of the first fit year (after burn-in)
#   fit_years          : integer vector of years (e.g., 2012:2019)
#   simulate_fit_years : function(par, y0, mdl, years) -> tibble(year, infections_total)
#   y_obs_df           : tibble with columns year, cases_reported_total
#   par_init_grid      : tibble of initial guesses (columns: beta, rho, k)
#   bounds             : list with lower/upper named numeric vectors for (beta, rho, k)
#   workers            : integer number of parallel workers
#
# Returns:
#   list(estimates = tibble(...), best = named vector, fit_series = tibble(...), details = list(...))
# -----------------------------------------------------------------------------
# R/fit_nb.R  (replace fit_nb_to_cases with this)
# -----------------------------------------------------------------------------
fit_nb_to_cases <- function(
    mdl,
    y0_burn,
    fit_years,
    simulate_fit_years,
    y_obs_df,
    par_init_grid,
    # Choose finite, conservative bounds so L-BFGS-B is stable
    bounds = list(
      lower = c(beta = 1e-2,  rho = 1e-8,  k = 1e-4),
      upper = c(beta = 5e3,   rho = 0.3,   k = 1e+5)
    ),
    workers = max(1L, parallel::detectCores() - 1L),
    sigma, gamma, VE  # pass the fixed epi params explicitly
) {
  stopifnot(all(c("year","cases_reported_total") %in% names(y_obs_df)))
  stopifnot(all(c("beta","rho","k") %in% names(par_init_grid)))
  stopifnot(all(c("beta","rho","k") %in% names(bounds$lower)))
  stopifnot(all(c("beta","rho","k") %in% names(bounds$upper)))
  
  # Align the data to the fit window (ascending years)
  y_obs_df <- y_obs_df %>% dplyr::filter(year %in% fit_years) %>% dplyr::arrange(year)
  years    <- as.integer(y_obs_df$year)
  y_obs    <- as.numeric(y_obs_df$cases_reported_total)
  
  # Parallel plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  
  # Objective wrapper for optim
  obj_fun <- function(theta) {
    beta <- theta[["beta"]]; rho <- theta[["rho"]]; k <- theta[["k"]]
    
    # Quick bound guard 
    if (any(!is.finite(c(beta, rho, k)))) return(Inf)
    if (beta < bounds$lower["beta"] || beta > bounds$upper["beta"]) return(Inf)
    if (rho  < bounds$lower["rho"]  || rho  > bounds$upper["rho"])  return(Inf)
    if (k    < bounds$lower["k"]    || k    > bounds$upper["k"])    return(Inf)
    
    # Simulate expected means for these parameters
    mu_hat <- try(
      simulate_mu_series(beta = beta, rho = rho,
                         y0_burn = y0_burn, mdl = mdl, fit_years = years,
                         simulate_fit_years = simulate_fit_years,
                         sigma = sigma, gamma = gamma, VE = VE),
      silent = TRUE
    )
    if (inherits(mu_hat, "try-error") || any(!is.finite(mu_hat))) return(Inf)
    if (length(mu_hat) != length(y_obs)) return(Inf)
    
    # NB negative log-likelihood
    nb_negloglik(y = y_obs, mu = mu_hat, k = k)
  }
  
  run_one_start <- function(start_row) {
    theta0 <- c(
      beta = as.numeric(start_row[["beta"]]),
      rho  = as.numeric(start_row[["rho"]]),
      k    = as.numeric(start_row[["k"]])
    )
    theta0 <- pmax(theta0, bounds$lower); theta0 <- pmin(theta0, bounds$upper)
    
    opt <- try(stats::optim(
      par     = theta0,
      fn      = obj_fun,
      method  = "L-BFGS-B",
      lower   = bounds$lower[c("beta","rho","k")],
      upper   = bounds$upper[c("beta","rho","k")],
      control = list(maxit = 200, factr = 1e7)
    ), silent = TRUE)
    
    if (inherits(opt, "try-error") ||
        any(!is.finite(unlist(opt$par))) ||
        !is.finite(opt$value)) {
      return(tibble::tibble(
        beta = theta0["beta"], rho = theta0["rho"], k = theta0["k"],
        nll = Inf, converged = FALSE, method = "L-BFGS-B"
      ))
    }
    
    tibble::tibble(
      beta = opt$par["beta"], rho = opt$par["rho"], k = opt$par["k"],
      nll  = as.numeric(opt$value),
      converged = isTRUE(opt$convergence == 0),
      method    = "L-BFGS-B"
    )
  }
  
  # Run all starts in parallel
  est_list <- furrr::future_map(
    .x = split(par_init_grid, seq_len(nrow(par_init_grid))),
    .f = run_one_start,
    .options = furrr::furrr_options(seed = TRUE)
  )
  # Combine all starts
  estimates <- dplyr::bind_rows(est_list)
  
  # Keep only finite rows (parameters and objective)
  estimates <- estimates %>%
    dplyr::filter(
      is.finite(beta), is.finite(rho), is.finite(k), is.finite(nll)
    ) %>%
    dplyr::arrange(nll, dplyr::desc(converged))
  
  if (nrow(estimates) == 0) {
    stop("All starts failed (no finite estimates). Consider widening bounds or changing starts.")
  }
  
  # Prefer converged; else take the smallest NLL nonetheless
  has_conv <- estimates %>% dplyr::filter(converged)
  best_row <- if (nrow(has_conv) > 0) has_conv[1, ] else estimates[1, ]
  
  best <- c(beta = as.numeric(best_row$beta),
            rho  = as.numeric(best_row$rho),
            k    = as.numeric(best_row$k))
  
  
  mu_best <- try(
    simulate_mu_series(
      beta = best["beta"], rho = best["rho"],
      y0_burn = y0_burn, mdl = mdl, fit_years = years,
      simulate_fit_years = simulate_fit_years,
      sigma = sigma, gamma = gamma, VE = VE
    ),
    silent = TRUE
  )
  if (inherits(mu_best, "try-error") || any(!is.finite(mu_best))) {
    stop("Simulation failed at the chosen best parameters. This suggests all starts were poor; ",
         "consider widening bounds or changing starting values.")
  }
  
  
  fit_series <- tibble::tibble(
    year   = years,
    obs    = y_obs,
    mu_hat = mu_best,
    resid  = y_obs - mu_best
  )
  
  list(
    estimates  = estimates,
    best       = best,
    fit_series = fit_series,
    details    = list(bounds = bounds, workers = workers, fit_years = years)
  )
}

# -----------------------------------------------------------------------------
# Helper: run the model to produce mu_t (expected reported counts) per year
# -----------------------------------------------------------------------------
# Inputs:
#   beta, rho         : numeric scalars to evaluate
#   y0_burn           : state vector at the start of the first fit year (already burned in)
#   mdl               : model object from build_model()
#   fit_years         : integer vector of years (ascending)
#   simulate_fit_years: function(par, y0, mdl, years) -> tibble(year, infections_total)
#   sigma, gamma, VE  : fixed epidemiological parameters (used by simulator)
# Returns:
#   numeric vector mu_t (same length as fit_years), where mu_t = rho * infections_total_t
simulate_mu_series <- function(beta, rho, y0_burn, mdl, fit_years,
                               simulate_fit_years, sigma, gamma, VE) {
  stopifnot(is.finite(beta), is.finite(rho), beta >= 0, rho >= 0)
  # Build parameter list used by the provided simulator
  pr_here <- list(sigma = sigma, gamma = gamma, beta = beta, rho = rho, VE = VE, eps = 1e-8)
  
  sim_df <- simulate_fit_years(par = pr_here, y0 = y0_burn, mdl = mdl, years = fit_years)
  stopifnot(all(c("year","infections_total") %in% names(sim_df)))
  
  sim_df <- sim_df %>%
    dplyr::filter(year %in% fit_years) %>%
    dplyr::arrange(year) %>%
    dplyr::mutate(mu = rho * infections_total)
  
  as.numeric(sim_df$mu)
}

# -----------------------------------------------------------------------------
# Helper: Negative Binomial negative log-likelihood
# -----------------------------------------------------------------------------
# Parameterization: NB(size = k, mu = mu_t)
#   Var(Y_t) = mu_t + mu_t^2 / k  with k > 0
# Inputs:
#   y  : non-negative integer vector (observed counts)
#   mu : non-negative numeric vector (expected means; same length as y)
#   k  : positive overdispersion (>0). Larger k -> closer to Poisson.
# Returns:
#   scalar negative log-likelihood (finite if inputs are valid)
nb_negloglik <- function(y, mu, k) {
  stopifnot(length(y) == length(mu))
  if (!is.finite(k) || k <= 0) return(Inf)
  
  # Guardrails: clip tiny/negative means to small positive to avoid log(0)
  mu <- pmax(mu, 1e-12)
  
  # dnbinom is vectorized and numerically stable for size/mu parameterization
  nll <- -sum(stats::dnbinom(x = y, size = k, mu = mu, log = TRUE))
  
  if (!is.finite(nll)) Inf else nll
}

