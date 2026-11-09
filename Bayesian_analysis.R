# ==============================================================================
# BAYESIAN & MONTE CARLO SIMULATION SCRIPT - MAIZE NITROGEN FIXATION STUDY
# ==============================================================================
# NOTE FOR REPRODUCIBILITY:
# This script relies on the 'maize_crops' dataset. Please note that the column 
# 'Valore_conc_azotata' (representing the estimated monetary cost of nitrogen 
# fertilization) is not natively present in the raw RICA microdata files. 
# As detailed in the README, this variable was pre-calculated 
# by processing data from a separate source file named 'RICA_Fertilizzanti.xlsx' 
# (Equations 1), whose microdata are subject to statistical confidentiality 
# and third-party access restrictions.

# ==============================================================================
# BAYESIAN & MONTE CARLO SIMULATION SCRIPT - MAIZE NITROGEN FIXATION STUDY
# ==============================================================================
# The script is divided in two parts: 
# 1. National Analysis (Aggregated across target regions)
# 2. Regional Analysis (Disaggregated by target regions and years)
# ==============================================================================

# Libraries used
library(dplyr)       # Data manipulation and creation of data frames
library(truncnorm)   # Truncated normal distributions
library(triangle)    # Triangular distribution functions (min, max, mean)
library(MASS)        # Advanced statistics, used here for mvrnorm (multivariate normal vectors)
library(coda)        # Monte Carlo simulation analysis


# ==============================================================================
# GLOBAL CONFIGURATION & COLUMN MAPPING SETUP
# ==============================================================================

# Define column references used consistently across the entire script
col_fertilizer <- "Valore_conc_azotata"
col_seeds <- "Sementi"
col_hectares <- "SUPERFICIE_UTIL"
col_region <- "Regione"
col_year <- "Anno"

# Define global simulation parameters and target correlation
set.seed(42)
n_sim <- 10000
indifference_threshold <- 30
target_final_cor <- 0.4

target_regions <- c(
  "Piemonte",
  "Emilia Romagna",
  "Veneto",
  "Friuli Venezia Giulia",
  "Lombardia"
)


# ==============================================================================
# PART 1: NATIONAL ANALYSIS 
# ==============================================================================

# Filter dataset for target regions
maize_crops_region_calibrated <- maize_crops %>%
  filter(!!sym(col_region) %in% target_regions)

available_years_national <- sort(unique(maize_crops_region_calibrated[[col_year]]))

# Calibrate latent correlation structure using copula parameters (National)
calibrate_latent_rho_national <- function(target_cor, n_sim_calib = 20000) {
  
  f <- function(rho_lat) {
    temp_sigma <- matrix(
      c(1, rho_lat,
        rho_lat, 1),
      nrow = 2
    )
    
    temp_z <- MASS::mvrnorm(n_sim_calib, mu = c(0, 0), Sigma = temp_sigma)
    
    temp_u1 <- pnorm(temp_z[, 1])
    temp_u2 <- pnorm(temp_z[, 2])
    
    temp_x1 <- truncnorm::qtruncnorm(
      temp_u1,
      a = 0,
      b = 1,
      mean = 0.30,
      sd = 0.05
    )
    
    temp_x2 <- triangle::qtriangle(
      temp_u2,
      a = 0.05,
      b = 0.16,
      c = 0.10
    )
    
    cor(temp_x1, temp_x2)
  }
  
  uniroot(function(r) f(r) - target_cor,
          interval = c(-0.99, 0.99))$root
}

latent_rho_star_national <- calibrate_latent_rho_national(target_final_cor)

sigma_calibrated <- matrix(
  c(1, latent_rho_star_national,
    latent_rho_star_national, 1),
  nrow = 2
)

global_results_calibrated <- data.frame()

# Loop through each available year for probabilistic simulation (National)
for (current_year in available_years_national) {
  
  subset_data_national <- maize_crops_region_calibrated %>%
    filter(!!sym(col_year) == current_year)
  
  if (nrow(subset_data_national) < 5) next
  
  # Aggregate total costs and area
  total_fertilizer_cost <- sum(subset_data_national[[col_fertilizer]], na.rm = TRUE)
  total_seed_cost <- sum(subset_data_national[[col_seeds]], na.rm = TRUE)
  total_hectares <- sum(subset_data_national[[col_hectares]], na.rm = TRUE)
  
  if (total_hectares <= 0) next
  
  # Calculate weighted means per hectare
  mean_fertilizer_cost_per_ha <- total_fertilizer_cost / total_hectares
  mean_seed_cost_per_ha <- total_seed_cost / total_hectares
  
  # Simulate correlated variables using multivariate normal and copula distributions
  sim_z <- MASS::mvrnorm(n_sim, mu = c(0, 0), Sigma = sigma_calibrated)
  
  sim_u1 <- pnorm(sim_z[, 1])
  sim_u2 <- pnorm(sim_z[, 2])
  
  sim_nitrogen_efficiency <- truncnorm::qtruncnorm(
    sim_u1,
    a = 0,
    b = 1,
    mean = 0.30,
    sd = 0.05
  )
  
  sim_seed_price_increase <- triangle::qtriangle(
    sim_u2,
    a = 0.05,
    b = 0.16,
    c = 0.10
  )
  
  nitrogen_savings <- mean_fertilizer_cost_per_ha * sim_nitrogen_efficiency
  extra_seed_cost <- mean_seed_cost_per_ha * sim_seed_price_increase
  
  net_benefit <- nitrogen_savings - extra_seed_cost
  
  temp_results_national <- data.frame(
    Year = as.factor(current_year),
    Benefit = net_benefit
  )
  
  global_results_calibrated <- rbind(
    global_results_calibrated,
    temp_results_national
  )
  
  cat(paste0(
    current_year,
    ": Mean EUR/ha = ", round(mean(net_benefit), 2),
    " | Risk < 0 = ",
    round(mean(net_benefit < 0) * 100, 2),
    "%\n"
  ))
}

# Compute final statistical summaries and credible intervals (National)
annual_statistics_calibrated <- global_results_calibrated %>%
  group_by(Year) %>%
  summarise(
    Mean_Benefit = mean(Benefit),
    Loss_Risk_Perc = mean(Benefit < 0) * 100,
    Insufficiency_Risk_Perc = mean(Benefit < indifference_threshold) * 100,
    Credible_Min_95 = quantile(Benefit, 0.025),
    Credible_Max_95 = quantile(Benefit, 0.975),
    .groups = "drop"
  )


# ==============================================================================
# PART 2: REGIONAL ANALYSIS 
# ==============================================================================

global_results_regional <- data.frame()

# Filter dataset for target regions (Regional)
filtered_data_regional <- maize_crops %>%
  filter(!!sym(col_region) %in% target_regions)

available_years_regional <- sort(unique(filtered_data_regional[[col_year]]))

# Correlation calibration function using copula parameters (Regional)
calibrate_latent_rho_regional <- function(target, n_sim_calib = 20000) {
  
  f <- function(rho_lat) {
    sigma_matrix_reg <- matrix(c(1, rho_lat, rho_lat, 1), 2)
    
    z_matrix_reg <- MASS::mvrnorm(n_sim_calib, mu = c(0, 0), Sigma = sigma_matrix_reg)
    
    u1_reg <- pnorm(z_matrix_reg[, 1])
    u2_reg <- pnorm(z_matrix_reg[, 2])
    
    x1_reg <- truncnorm::qtruncnorm(u1_reg, a = 0, b = 1, mean = 0.30, sd = 0.05)
    x2_reg <- triangle::qtriangle(u2_reg, a = 0.05, b = 0.16, c = 0.10)
    
    cor(x1_reg, x2_reg)
  }
  
  uniroot(function(r) f(r) - target, interval = c(-0.99, 0.99))$root
}

# Calibrate latent correlation parameter (Regional)
latent_rho_star_regional <- calibrate_latent_rho_regional(target_final_cor, n_sim_calib = 20000)

sigma_matrix_regional <- matrix(
  c(1, latent_rho_star_regional,
    latent_rho_star_regional, 1),
  2
)

# MAIN LOOP: Iterate over regions and years for probabilistic simulation (Regional)
for (current_region in target_regions) {
  
  for (current_year in available_years_regional) {
    
    subset_data_regional <- filtered_data_regional %>%
      filter(!!sym(col_region) == current_region,
             !!sym(col_year) == current_year)
    
    if (nrow(subset_data_regional) < 5) next
    
    # Aggregate total costs and area
    total_fertilizer_cost_reg <- sum(subset_data_regional[[col_fertilizer]], na.rm = TRUE)
    total_seed_cost_reg <- sum(subset_data_regional[[col_seeds]], na.rm = TRUE)
    total_hectares_reg <- sum(subset_data_regional[[col_hectares]], na.rm = TRUE)
    
    if (total_hectares_reg <= 0) next
    
    # Calculate weighted means per hectare
    mean_fertilizer_cost_per_ha_reg <- total_fertilizer_cost_reg / total_hectares_reg
    mean_seed_cost_per_ha_reg <- total_seed_cost_reg / total_hectares_reg
    
    # Simulate correlated variables using multivariate normal and copula distributions
    z_matrix_reg <- MASS::mvrnorm(n_sim, mu = c(0, 0), Sigma = sigma_matrix_regional)
    
    u1_reg <- pnorm(z_matrix_reg[, 1])
    u2_reg <- pnorm(z_matrix_reg[, 2])
    
    sim_nitrogen_efficiency_reg <- truncnorm::qtruncnorm(
      u1_reg,
      a = 0,
      b = 1,
      mean = 0.30,
      sd = 0.05
    )
    
    sim_seed_price_increase_reg <- triangle::qtriangle(
      u2_reg,
      a = 0.05,
      b = 0.16,
      c = 0.10
    )
    
    # Calculate economic benefits and net returns
    nitrogen_savings_reg <- mean_fertilizer_cost_per_ha_reg * sim_nitrogen_efficiency_reg
    extra_seed_cost_reg <- mean_seed_cost_per_ha_reg * sim_seed_price_increase_reg
    
    net_benefit_reg <- nitrogen_savings_reg - extra_seed_cost_reg
    
    temp_results_regional <- data.frame(
      Region = current_region,
      Year = current_year,
      Benefit = net_benefit_reg
    )
    
    global_results_regional <- rbind(
      global_results_regional,
      temp_results_regional
    )
    
    cat(paste0(
      current_region, " - ", current_year,
      ": Mean EUR/ha = ", round(mean(net_benefit_reg), 2),
      " | Risk < 0 = ",
      round(mean(net_benefit_reg < 0) * 100, 2),
      "%\n"
    ))
    
  }
}

# Compute final regional and annual statistical summaries and credible intervals
regional_statistics <- global_results_regional %>%
  group_by(Region, Year) %>%
  summarise(
    Mean_Benefit = mean(Benefit),
    Loss_Risk_Perc = mean(Benefit < 0) * 100,
    Insufficiency_Risk_Perc = mean(Benefit < indifference_threshold) * 100,
    Credible_Min_95 = quantile(Benefit, 0.025),
    Credible_Max_95 = quantile(Benefit, 0.975),
    .groups = "drop"
  )
