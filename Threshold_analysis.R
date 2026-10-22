#This script aggregates farm-level simulation results to evaluate the macro-sectoral economic impact of adopting
#nitrogen-fixing maize. It structures the data across simulation scenarios and applies an adoption rule based on a
#minimum profitability threshold per hectare (set to 30 €/ha, with alternative analyses at 10 €/ha and 20 €/ha 
#obtained by modifying the adoption condition). By filtering for farms that meet this threshold, the script 
#computes adopting shares, total benefits, and the percentage increase in gross margin across the entire RICA 
#sample.
library(dplyr)
library(tidyr)

# 1. Convert to long format and zero out benefits below the threshold (< 30 €/ha)
maize_long <- maize_crops %>%
  pivot_longer(
    cols = c("Benefit_N_fix_00_Ha", "Benefit_N_fix_05_Ha", 
             "Benefit_N_fix_10_Ha", "Benefit_N_fix_16_Ha"),
    names_to = "Simulation",
    values_to = "Benefit_per_ha"
  ) %>%
  mutate(
    Simulation = case_when(
      grepl("00", Simulation) ~ "+0%",
      grepl("05", Simulation) ~ "+5%",
      grepl("10", Simulation) ~ "+10%",
      grepl("16", Simulation) ~ "+16%"
    ),
    # Adoption flag (1 if the 30 €/ha threshold is met or exceeded, 0 otherwise)
    Adopter = ifelse(Benefit_per_ha >= 30, 1, 0),
    
    # If the farm does not exceed the threshold, the benefit is set to 0 €
    Effective_Benefit_per_ha = ifelse(Adopter == 1, Benefit_per_ha, 0),
    Effective_Total_Benefit = Effective_Benefit_per_ha * SUPERFICIE_UTIL
  )

# 2. Calculate annual totals for the ENTIRE sample (fixed denominator)
annual_sample_totals <- maize_crops %>%
  group_by(Anno) %>%
  summarise(
    Total_Sample_Observations = n(),
    Total_Sample_Hectares = sum(SUPERFICIE_UTIL, na.rm = TRUE),
    
    # Total Gross Margin of ALL sample farms for the reference year (€)
    Total_Sample_Gross_Margin = sum(Gross_margin, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Aggregate by Year and Simulation with overall sectoral impact calculation
pivot_table <- maize_long %>%
  group_by(Anno, Simulation) %>%
  summarise(
    Adopting_Observations = sum(Adopter == 1, na.rm = TRUE),
    Adopting_Hectares = sum(SUPERFICIE_UTIL[Adopter == 1], na.rm = TRUE),
    Total_Benefit_Euro = sum(Effective_Total_Benefit, na.rm = TRUE),
    
    # Sum of total Gross Margin for adopting farms only (€)
    Adopting_Gross_Margin_Euro = sum(Gross_margin[Adopter == 1], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(annual_sample_totals, by = "Anno") %>%
  mutate(
    Pct_Adopting_Farms = (Adopting_Observations / Total_Sample_Observations) * 100,
    Pct_Adopting_Hectares = (Adopting_Hectares / Total_Sample_Hectares) * 100,
    
    # Mean benefit per hectare for adopting farms only (€/ha)
    Mean_Benefit_Adopters_per_ha = ifelse(Adopting_Hectares > 0, Total_Benefit_Euro / Adopting_Hectares, 0),
    
    # Percentage increase of Gross Margin calculated on the ENTIRE SAMPLE (macro-sectoral impact)
    Pct_Sectoral_Gross_Margin_Increase = (Total_Benefit_Euro / Total_Sample_Gross_Margin) * 100
  ) %>%
  select(
    Anno,
    Simulation,
    Adopting_Observations,
    Pct_Adopting_Farms,
    Adopting_Hectares,
    Pct_Adopting_Hectares,
    Total_Benefit_Euro,
    Mean_Benefit_Adopters_per_ha,
    Adopting_Gross_Margin_Euro,
    Total_Sample_Gross_Margin,
    Pct_Sectoral_Gross_Margin_Increase
  )