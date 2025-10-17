## ----setup, include=FALSE, warning=FALSE, message=FALSE----
# ---- setup ----
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
library(readr); library(dplyr); library(tidyr); library(stringr); library(deSolve); library(tibble)
set.seed(42)
verbose <- TRUE


## ----------------------------------------------
# ---- build_model ---------------------------------------------------------
# Creates a self-contained model object ("mdl") that holds:
# - structure (ages, sexes, sizes, index map)
# - lookups (mortality-by-year, pregnancy prevalence)
# - data tables (births, coverage)
# - helpers (cell_index)
#
# Inputs:
#   ages, sexes: vectors (e.g., 0:100, c("F","M"))
#   pop_sy: population at start_year by (sex, age), cols: year, sex, age, pop
#   mu_tbl: mortality hazards, cols: year, sex, age, mu (for all sim years)
#   births_df: cols year, births_total
#   coverage_df: cols year, routine_coverage, catchup_coverage
#   preg_prev: cols age (15..44), prev
build_model <- function(ages, sexes, pop_sy, mu_tbl, births_df, coverage_df, preg_prev) {

  # grid & sizes
  state_key <- expand.grid(sex = sexes, age = ages) |>
    dplyr::arrange(sex, age) |>
    dplyr::mutate(pos = dplyr::row_number())
  n_age <- length(ages); n_sex <- length(sexes); n_cells <- n_age * n_sex

  # index map (flat blocks)
  idx <- list(
    S = 1:n_cells,
    E = (n_cells + 1):(2*n_cells),
    I = (2*n_cells + 1):(3*n_cells),
    R = (3*n_cells + 1):(4*n_cells),
    C_F  = (4*n_cells + 1):(4*n_cells + n_age),
    C_tot = 4*n_cells + n_age + 1
  )

  # mortality lookup: year -> mu_vec ordered (sex-major, age-minor)
  build_mu_lookup <- function(mu_table) {
    years_avail <- sort(unique(mu_table$year))
    out <- setNames(vector("list", length(years_avail)), as.character(years_avail))
    for (yy in years_avail) {
      v <- mu_table |>
        dplyr::filter(year == yy) |>
        dplyr::arrange(sex, age) |>
        dplyr::pull(mu)
      stopifnot(length(v) == n_cells)
      out[[as.character(yy)]] <- v
    }
    out
  }
  mu_lookup <- build_mu_lookup(mu_tbl)

  # pregnancy prevalence vector aligned to ages (0..100)
  preg_prev_vec <- { v <- numeric(n_age); v[match(preg_prev$age, ages)] <- preg_prev$prev; v }

  # helper: linear cell index within a block
  cell_index <- function(sex, age) {
    s_idx <- match(sex, sexes) - 1L
    a_idx <- match(age, ages)  - 1L
    s_idx * n_age + (a_idx + 1L)
  }

  # which rows in a compartment block are female ages (first n_age)
  female_block_indices <- 1:n_age

  # starter population vector (sex-major, age-minor)
  pop0_vec <- pop_sy |>
    dplyr::arrange(sex, age) |>
    dplyr::pull(pop)

  list(
    # structure
    ages = ages, sexes = sexes, n_age = n_age, n_sex = n_sex, n_cells = n_cells,
    state_key = state_key, idx = idx,
    # lookups / data
    mu_lookup = mu_lookup,
    births_df = births_df,
    coverage_df = coverage_df,
    preg_prev_vec = preg_prev_vec,
    # helpers
    cell_index = cell_index,
    female_block_indices = female_block_indices,
    # initial pop vector (for building y0)
    pop0_vec = pop0_vec
  )
}


## ----------------------------------------------
# ---- build_initial_state -------------------------------------------------
# Builds y0 using mdl’s structure. Seeds 1 infection per cell (if pop>0).
build_initial_state <- function(mdl, seed_I_per_cell = 1L) {
  with(mdl, {
    I0 <- ifelse(pop0_vec > 0, pmin(seed_I_per_cell, pop0_vec), 0)
    E0 <- rep(0, n_cells); R0 <- rep(0, n_cells)
    S0 <- pmax(pop0_vec - I0, 0)
    y0 <- c(S = S0, E = E0, I = I0, R = R0, C_F = rep(0, n_age), C_tot = 0)
    y0
  })
}


## ----------------------------------------------
# ---- seir_ode_mdl -------------------------------------------------------
# ODE uses:
#   - year-specific mu via pr$mu_vec
#   - structure & indices via pr$mdl
seir_ode_mdl <- function(t, y, pr) {
  if (is.null(pr$mdl)) stop("seir_ode_mdl: parms$mdl is NULL")
  if (is.null(pr$mu_vec)) stop("seir_ode_mdl: parms$mu_vec is NULL")
  mdl <- pr$mdl
  
  with(mdl, {
    S <- y[idx$S]; E <- y[idx$E]; I <- y[idx$I]; R <- y[idx$R]
    
    # Hard guards: fail fast if any non-finite values appear
    if (!all(is.finite(S))) stop("seir_ode_mdl: non-finite values in S")
    if (!all(is.finite(E))) stop("seir_ode_mdl: non-finite values in E")
    if (!all(is.finite(I))) stop("seir_ode_mdl: non-finite values in I")
    if (!all(is.finite(R))) stop("seir_ode_mdl: non-finite values in R")
    
    N_tot <- sum(S + E + I + R)
    I_tot <- sum(I)
    if (!is.finite(N_tot)) stop("seir_ode_mdl: N_tot is non-finite")
    
    # Tiny importation floor allowed via pr$eps (default 0)
    eps <- if (!is.null(pr$eps) && is.finite(pr$eps)) pr$eps else 0
    lam_core <- if (N_tot > 0) pr$beta * (I_tot / N_tot) else 0
    if (!is.finite(lam_core)) stop("seir_ode_mdl: lambda core is non-finite")
    lambda <- lam_core + eps
    
    dS <- -lambda * S
    dE <-  lambda * S - pr$sigma * E
    dI <-  pr$sigma * E - pr$gamma * I
    dR <-  pr$gamma * I
    
    mu_here <- pr$mu_vec
    if (length(mu_here) != length(S)) {
      stop(sprintf("seir_ode_mdl: mu_vec length mismatch: got %d, need %d",
                   length(mu_here), length(S)))
    }
    if (!all(is.finite(mu_here))) stop("seir_ode_mdl: non-finite mortality hazard")
    
    dS <- dS - mu_here * S
    dE <- dE - mu_here * E
    dI <- dI - mu_here * I
    dR <- dR - mu_here * R
    
    S_F    <- S[female_block_indices]
    dC_F   <- lambda * S_F
    dC_tot <- lambda * sum(S)
    
    list(c(dS, dE, dI, dR, dC_F, dC_tot))
  })
}



## ----------------------------------------------
# ---- events_mdl ----------------------------------------------------------
apply_births_mdl <- function(y, yr, mdl, births_tbl, share_F, share_M) {
  B <- births_tbl |> dplyr::filter(year == yr) |> dplyr::pull(births_total)
  if (length(B) == 0 || is.na(B) || B <= 0) return(y)
  iF0 <- mdl$idx$S[ mdl$cell_index("F", 0L) ]
  iM0 <- mdl$idx$S[ mdl$cell_index("M", 0L) ]
  y[iF0] <- y[iF0] + B * share_F
  y[iM0] <- y[iM0] + B * share_M
  y
}

apply_ageing_mdl <- function(y, mdl) {
  shift_one_block <- function(block) {
    new_block <- numeric(length(block))
    for (sx_i in seq_along(mdl$sexes)) {
      base <- (sx_i - 1L) * mdl$n_age
      from <- base + 1L:(mdl$n_age - 1L)
      to   <- base + 2L:mdl$n_age
      new_block[to] <- new_block[to] + block[from]
      new_block[base + mdl$n_age] <- new_block[base + mdl$n_age] + block[base + mdl$n_age]
    }
    new_block
  }
  y[mdl$idx$S] <- shift_one_block(y[mdl$idx$S])
  y[mdl$idx$E] <- shift_one_block(y[mdl$idx$E])
  y[mdl$idx$I] <- shift_one_block(y[mdl$idx$I])
  y[mdl$idx$R] <- shift_one_block(y[mdl$idx$R])
  y
}

apply_routine_mdl <- function(y, yr, mdl, cov_tbl, VE) {
  cov <- cov_tbl |> dplyr::filter(year == yr) |> dplyr::pull(routine_coverage)
  if (length(cov) == 0 || is.na(cov) || cov <= 0) return(y)
  eff <- max(0, min(cov * VE, 1))
  for (sx in mdl$sexes) {
    iS <- mdl$idx$S[ mdl$cell_index(sx, 1L) ]
    iR <- mdl$idx$R[ mdl$cell_index(sx, 1L) ]
    move <- eff * y[iS]
    y[iS] <- y[iS] - move; y[iR] <- y[iR] + move
  }
  y
}

apply_catchup_mdl <- function(y, yr, mdl, cov_tbl, VE) {
  cov <- cov_tbl |> dplyr::filter(year == yr) |> dplyr::pull(catchup_coverage)
  if (length(cov) == 0 || is.na(cov) || cov <= 0) return(y)
  eff <- max(0, min(cov * VE, 1))
  for (sx in mdl$sexes) {
    cells <- vapply(1:14, function(a) mdl$cell_index(sx, a), integer(1))
    iS <- mdl$idx$S[cells]; iR <- mdl$idx$R[cells]
    move <- eff * y[iS]
    y[iS] <- y[iS] - move; y[iR] <- y[iR] + move
  }
  y
}


## ----------------------------------------------
# ---- simulate_one_year_mdl ----------------------------------------------
simulate_one_year_mdl <- function(y_start, year_int, pr, mdl, share_F, share_M) {
  # add year-specific mortality & the model object into parms
  mu_year <- mdl$mu_lookup[[as.character(year_int)]]
  if (is.null(mu_year)) stop("No mortality vector for year ", year_int)
  pr_year <- pr; pr_year$mu_vec <- mu_year; pr_year$mdl <- mdl

  # integrate one calendar year
  times <- seq(year_int, year_int + 1, by = 1/365)
  sol <- deSolve::ode(y = y_start, times = times, func = seir_ode_mdl, parms = pr_year, method = "rk4")

  # pre-events end-of-year state
  y_end <- as.numeric(sol[nrow(sol), -1])

  # events
  y_ev <- y_end
  y_ev <- apply_ageing_mdl  (y_ev, mdl)
  y_ev <- apply_births_mdl (y_ev, year_int, mdl, mdl$births_df, share_F, share_M)
  y_ev <- apply_routine_mdl (y_ev, year_int, mdl, mdl$coverage_df, pr$VE)
  y_ev <- apply_catchup_mdl (y_ev, year_int, mdl, mdl$coverage_df, pr$VE)

  list(trajectory = sol, y_next = y_ev)
}


## ----------------------------------------------
# ---- crs_from_infections_mdl --------------------------------------------
crs_from_infections_mdl <- function(inf_F_by_age, mdl, r, w) {
  stopifnot(length(inf_F_by_age) == mdl$n_age)
  r_eff <- r$r1 * w["w1"] + r$r2 * w["w2"] + r$r3 * w["w3"]
  sum(inf_F_by_age * mdl$preg_prev_vec * r_eff)
}


## ----------------------------------------------
# ---- run_horizon_mdl -----------------------------------------------------
run_horizon_mdl <- function(y_start, start_year, n_years, pr, mdl, share_F, share_M, crs_parms, tri_weights) {
  yrs <- seq.int(from = start_year, length.out = n_years)
  out <- tibble::tibble(
    year = integer(),
    infections_total = double(),
    reported_cases   = double(),
    crs              = double()
  )
  y_curr <- y_start
  for (yr in yrs) {
    # reset counters at year start
    y_curr[mdl$idx$C_F]  <- 0
    y_curr[mdl$idx$C_tot] <- 0

    res   <- simulate_one_year_mdl(y_curr, yr, pr, mdl, share_F, share_M)
    y_pre <- as.numeric(res$trajectory[nrow(res$trajectory), -1])

    C_F_end   <- y_pre[mdl$idx$C_F]
    C_tot_end <- as.numeric(y_pre[mdl$idx$C_tot])

    out <- dplyr::bind_rows(out, tibble::tibble(
      year = yr,
      infections_total = C_tot_end,
      reported_cases   = pr$rho * C_tot_end,
      crs              = crs_from_infections_mdl(C_F_end, mdl, r = crs_parms, w = tri_weights)
    ))

    y_curr <- res$y_next
  }
  list(results = out, y_final = y_curr)
}


## ----------------------------------------------
# ---- burn_in_mdl ---------------------------------------------------------
make_stationary_births <- function(births_base_tbl, years, ref_year) {
  B_ref <- births_base_tbl |> dplyr::filter(year == ref_year) |> dplyr::pull(births_total)
  stopifnot(length(B_ref) == 1, is.finite(B_ref))
  tibble::tibble(year = years, births_total = rep(B_ref, length(years)))
}

make_zero_coverage <- function(years) {
  tibble::tibble(year = years, routine_coverage = 0, catchup_coverage = NA_real_)
}

# Replacement-births variant: births = continuous deaths to keep N ~ constant
simulate_one_year_custom_mdl <- function(
    y_start, year_int, pr, mdl,
    coverage_tbl, births_tbl,
    share_F, share_M,
    mu_year_override = NULL,
    replacement_births = FALSE      # <— new flag
) {
  mu_year <- if (!is.null(mu_year_override)) mu_year_override else mdl$mu_lookup[[as.character(year_int)]]
  if (is.null(mu_year)) stop("No mortality vector for year ", year_int)
  pr_year <- pr; pr_year$mu_vec <- mu_year; pr_year$mdl <- mdl
  
  # integrate one calendar year
  times <- seq(year_int, year_int + 1, by = 1/365)
  sol <- deSolve::ode(y = y_start, times = times, func = seir_ode_mdl, parms = pr_year, method = "rk4")
  
  # end-of-year state (pre events)
  y_end <- as.numeric(sol[nrow(sol), -1])
  
  # compute births to add
  if (replacement_births) {
    N_start <- sum(y_start[mdl$idx$S] + y_start[mdl$idx$E] + y_start[mdl$idx$I] + y_start[mdl$idx$R])
    N_end   <- sum(y_end  [mdl$idx$S] + y_end  [mdl$idx$E] + y_end  [mdl$idx$I] + y_end  [mdl$idx$R])
    B_here  <- max(0, N_start - N_end)  # replace continuous deaths; never negative
    births_tbl_here <- tibble::tibble(year = year_int, births_total = B_here)
  } else {
    births_tbl_here <- births_tbl
  }
  
  # apply discrete events
  y_ev <- y_end
  y_ev <- apply_ageing_mdl  (y_ev, mdl)
  y_ev <- apply_births_mdl (y_ev, year_int, mdl, births_tbl_here, share_F, share_M)
  y_ev <- apply_routine_mdl (y_ev, year_int, mdl, coverage_tbl, pr$VE)
  y_ev <- apply_catchup_mdl (y_ev, year_int, mdl, coverage_tbl, pr$VE)
  
  list(trajectory = sol, y_next = y_ev)
}


burn_in_mdl <- function(y_start, start_year, K, pr, mdl, ref_year, share_F, share_M) {
  if (K <= 0) return(y_start)
  yrs_burn    <- seq.int(from = start_year - K, length.out = K)
  cov_zero    <- make_zero_coverage(yrs_burn)
  births_stat <- make_stationary_births(mdl$births_df, yrs_burn, ref_year = ref_year)
  mu_ref      <- mdl$mu_lookup[[as.character(ref_year)]]
  if (is.null(mu_ref)) stop("No mortality for ref_year ", ref_year)
  
  y <- y_start
  for (yr in yrs_burn) {
    y[mdl$idx$C_F] <- 0; y[mdl$idx$C_tot] <- 0
    res <- simulate_one_year_custom_mdl(
      y_start  = y, year_int = yr, pr = pr, mdl = mdl,
      coverage_tbl = cov_zero,
      births_tbl   = births_stat,
      share_F = share_F, share_M = share_M,
      mu_year_override = mu_ref,
      replacement_births = TRUE     # <— keep N ~ constant during burn-in
    )
    y <- res$y_next
  }
  y
}

# ----------------------------------------------
# ---- run_stationary_novax -----------------------------------------------
# Simulate a set of calendar years with:
#   - stationary demography (mortality fixed at ref_year),
#   - replacement births to keep N approximately stable,
#   - zero vaccination (routine=0, no catch-up).
#
# Inputs:
#   y_start   : state vector at the start of the first year in `years`
#   years     : integer vector of years to simulate (e.g., 2012:2019)
#   pr        : parameter list (sigma, gamma, beta, rho, VE, eps)
#   mdl       : model object from build_model()
#   ref_year  : year whose mortality is used for all steps (stationary demography)
#   share_F/M : sex ratio splitting for births
#
# Returns:
#   list(results = tibble with yearly infections_total, reported_cases, crs,
#                 plus optional diagnostics if desired),
#        y_final = state vector at end of the last simulated year)
run_stationary_novax <- function(y_start, years, pr, mdl, ref_year, share_F, share_M) {
  stopifnot(length(years) > 0, is.numeric(years))
  yrs <- sort(unique(as.integer(years)))
  
  # Zero-coverage table for the provided years
  cov_zero <- tibble::tibble(
    year = yrs,
    routine_coverage = 0,
    catchup_coverage = NA_real_
  )
  
  # Use stationary births each year (constant on the ref year’s total),
  # then override with "replacement births" inside the one-year simulator.
  births_stat <- make_stationary_births(mdl$births_df, yrs, ref_year = ref_year)
  
  # Fixed mortality from ref_year for all years
  mu_ref <- mdl$mu_lookup[[as.character(ref_year)]]
  if (is.null(mu_ref)) stop("run_stationary_novax: no mortality for ref_year = ", ref_year)
  
  out <- tibble::tibble(
    year = integer(),
    infections_total = double(),
    reported_cases   = double(),
    crs              = double()
  )
  
  y_curr <- y_start
  
  for (yr in yrs) {
    # Reset annual accumulators before integrating
    y_curr[mdl$idx$C_F]   <- 0
    y_curr[mdl$idx$C_tot] <- 0
    
    # Integrate one calendar year with:
    #  - fixed mortality (mu_year_override = mu_ref),
    #  - zero vaccination (cov_zero),
    #  - replacement_births = TRUE to keep N ~ stable
    res <- simulate_one_year_custom_mdl(
      y_start  = y_curr,
      year_int = yr,
      pr       = pr,
      mdl      = mdl,
      coverage_tbl = cov_zero[cov_zero$year == yr, , drop = FALSE],
      births_tbl   = births_stat[births_stat$year == yr, , drop = FALSE],
      share_F = share_F, share_M = share_M,
      mu_year_override   = mu_ref,
      replacement_births = TRUE
    )
    
    # Read accumulators BEFORE events are applied (already built into simulate_one_year_custom_mdl flow)
    y_pre <- as.numeric(res$trajectory[nrow(res$trajectory), -1])
    C_F_end   <- y_pre[mdl$idx$C_F]
    C_tot_end <- as.numeric(y_pre[mdl$idx$C_tot])
    
    # Compute CRS with existing helper
    crs_year <- crs_from_infections_mdl(C_F_end, mdl, r = pr$crs_parms %||% list(r1 = 0.90, r2 = 0.25, r3 = 0.01),
                                        w = pr$tri_weights %||% c(w1 = 1/3, w2 = 1/3, w3 = 1/3))
    
    out <- dplyr::bind_rows(out, tibble::tibble(
      year = yr,
      infections_total = C_tot_end,
      reported_cases   = if (!is.null(pr$rho) && is.finite(pr$rho)) pr$rho * C_tot_end else NA_real_,
      crs              = crs_year
    ))
    
    # Move to next year (post-events state)
    y_curr <- res$y_next
  }
  
  list(results = out, y_final = y_curr)
}

