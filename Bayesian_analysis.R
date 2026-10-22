# The following script performs the Bayesian analysis developed for this study. 
# The script is divided in two parts: The first part runs an analysis starting from 
# the "maize_crops" data frame (the file containing the RICA microdata), which computes 
# probabilities in an aggregated manner. 
# The second part repeats the same analysis, but breaks down each result by the considered regions. 
# Please ensure that R has the columns related to years and farm ID codes set to numeric format.

#library used
library(dplyr)        # Data manipulation and creation of data frames
library(truncnorm)    # Truncated normal distributions
library(triangle)     # Triangular distribution functions (min, max, mean)
library(MASS)         # Advanced statistics, used here for mvrnorm (multivariate normal vectors)
library(coda)         # Monte Carlo simulation analysis

#National Analysis ----

set.seed(42)

# Define simulation parameters and target correlation
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

# Filter dataset for target regions
maize_crops_region_calibrated <- maize_crops %>%
  filter(Region %in% target_regions)

available_years <- sort(unique(maize_crops_region_calibrated$Year))

# Calibrate latent correlation structure using copula parameters
calibrate_latent_rho <- function(target_cor, n_sim = 20000) {
  
  f <- function(rho_lat) {
    temp_sigma <- matrix(
      c(1, rho_lat,
        rho_lat, 1),
      nrow = 2
    )
    
    temp_z <- MASS::mvrnorm(n_sim, mu = c(0, 0), Sigma = temp_sigma)
    
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

latent_rho_star <- calibrate_latent_rho(target_final_cor)

sigma_calibrated <- matrix(
  c(1, latent_rho_star,
    latent_rho_star, 1),
  nrow = 2
)

global_results_calibrated <- data.frame()

# Loop through each available year for probabilistic simulation
for (current_year in available_years) {
  
  subset_data <- maize_crops_region_calibrated %>%
    filter(Year == current_year)
  
  if (nrow(subset_data) < 5) next
  
  # Aggregate total costs and area
  total_fertilizer_cost <- sum(subset_data[[col_fertilizer]], na.rm = TRUE)
  total_seed_cost <- sum(subset_data[[col_seeds]], na.rm = TRUE)
  total_hectares <- sum(subset_data[[col_hectares]], na.rm = TRUE)
  
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
  
  temp_results <- data.frame(
    Year = as.factor(current_year),
    Benefit = net_benefit
  )
  
  global_results_calibrated <- rbind(
    global_results_calibrated,
    temp_results
  )
  
  cat(paste0(
    current_year,
    ": Mean EUR/ha = ", round(mean(net_benefit), 2),
    " | Risk < 0 = ",
    round(mean(net_benefit < 0) * 100, 2),
    "%\n"
  ))
}

# Compute final statistical summaries and credible intervals
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

####################################################################################
##################################################################################

# REGIONAL ANALYSIS SCRIPT

## SIMULATION PARAMETERS
set.seed(42)
n_sim <- 10000
indifference_threshold <- 30
target_cor <- 0.4

global_results_regional <- data.frame()

# COLUMN CONFIGURATION
col_fertilizer <- "Valore_conc_azotata"
col_seeds <- "Sementi"
col_hectares <- "SUPERFICIE_UTIL"
col_region <- "Regione"

target_regions <- c(
  "Piemonte",
  "Lombardia",
  "Veneto",
  "Friuli Venezia Giulia",
  "Emilia Romagna"
)

# Filter dataset for target regions
filtered_data <- maize_crops %>%
  filter(!!sym(col_region) %in% target_regions)

available_years <- sort(unique(filtered_data$Anno))

# Correlation calibration function using copula parameters
calibrate_latent_rho <- function(target, n_sim = 20000) {
  
  f <- function(rho_lat) {
    sigma_matrix <- matrix(c(1, rho_lat, rho_lat, 1), 2)
    
    z_matrix <- MASS::mvrnorm(n_sim, mu = c(0, 0), Sigma = sigma_matrix)
    
    u1 <- pnorm(z_matrix[, 1])
    u2 <- pnorm(z_matrix[, 2])
    
    x1 <- truncnorm::qtruncnorm(u1, a = 0, b = 1, mean = 0.30, sd = 0.05)
    x2 <- triangle::qtriangle(u2, a = 0.05, b = 0.16, c = 0.10)
    
    cor(x1, x2)
  }
  
  uniroot(function(r) f(r) - target, interval = c(-0.99, 0.99))$root
}

# Calibrate latent correlation parameter
latent_rho_star <- calibrate_latent_rho(target_cor, n_sim = 20000)

sigma_matrix <- matrix(
  c(1, latent_rho_star,
    latent_rho_star, 1),
  2
)

# MAIN LOOP: Iterate over regions and years for probabilistic simulation
for (current_region in target_regions) {
  
  for (current_year in available_years) {
    
    subset_data <- filtered_data %>%
      filter(Regione == current_region,
             Anno == current_year)
    
    if (nrow(subset_data) < 5) next
    
    # Aggregate total costs and area
    total_fertilizer_cost <- sum(subset_data[[col_fertilizer]], na.rm = TRUE)
    total_seed_cost <- sum(subset_data[[col_seeds]], na.rm = TRUE)
    total_hectares <- sum(subset_data[[col_hectares]], na.rm = TRUE)
    
    if (total_hectares <= 0) next
    
    # Calculate weighted means per hectare
    mean_fertilizer_cost_per_ha <- total_fertilizer_cost / total_hectares
    mean_seed_cost_per_ha <- total_seed_cost / total_hectares
    
    # Simulate correlated variables using multivariate normal and copula distributions
    z_matrix <- MASS::mvrnorm(n_sim, mu = c(0, 0), Sigma = sigma_matrix)
    
    u1 <- pnorm(z_matrix[, 1])
    u2 <- pnorm(z_matrix[, 2])
    
    sim_nitrogen_efficiency <- truncnorm::qtruncnorm(
      u1,
      a = 0,
      b = 1,
      mean = 0.30,
      sd = 0.05
    )
    
    sim_seed_price_increase <- triangle::qtriangle(
      u2,
      a = 0.05,
      b = 0.16,
      c = 0.10
    )
    
    # Calculate economic benefits and net returns
    nitrogen_savings <- mean_fertilizer_cost_per_ha * sim_nitrogen_efficiency
    extra_seed_cost <- mean_seed_cost_per_ha * sim_seed_price_increase
    
    net_benefit <- nitrogen_savings - extra_seed_cost
    
    temp_results <- data.frame(
      Region = current_region,
      Year = current_year,
      Benefit = net_benefit
    )
    
    global_results_regional <- rbind(
      global_results_regional,
      temp_results
    )
    
    cat(paste0(
      current_region, " - ", current_year,
      ": Mean EUR/ha = ", round(mean(net_benefit), 2),
      " | Risk < 0 = ",
      round(mean(net_benefit < 0) * 100, 2),
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