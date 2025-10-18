# ---- config.R -------------------------------------------------------------
# Central knobs used by analysis scripts.

# ---- config.R (policy horizon) ----
set.seed(42)
verbose <- TRUE

start_year    <- 2025
horizon_years <- 25  # <-- set 10 for short-term; 20–25 to see routine's cohort effect

sim_years     <- start_year:(start_year + horizon_years - 1)

# Demographic structure
ages  <- 0:100
sexes <- c("F","M")

# Epidemiology (placeholders; beta, rho will be fitted later)
latent_period_years     <- 14/365
infectious_period_years <-  7/365
VE   <- 0.95
beta <- 7.5*(1/infectious_period_years)   # from literature
rho  <- 0.01   # from literature

# CRS risks and trimester weights (can vary in sensitivity)
crs_parms   <- list(r1 = 0.90, r2 = 0.25, r3 = 0.01)
tri_weights <- c(w1 = 1/3, w2 = 1/3, w3 = 1/3)

# Sex ratio at birth (≈1.05 → ~51.2% male)
SRB     <- 1.05
share_F <- 1 / (1 + SRB)
share_M <- SRB / (1 + SRB)
