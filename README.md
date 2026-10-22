# nitrogen-fixing-maize-simulation-italy
This repository hosts the R source code and analytical scripts developed to evaluate the economic feasibility and farm-level adoption dynamics of introducing a novel nitrogen-fixing maize variety (analogous to the Sierra Mixe trait) within Italian agricultural systems. The computational framework integrates traditional partial budget modeling with advanced probabilistic and Bayesian techniques to rigorously account for agronomic and economic uncertainties.

Bayesian_Analysis.R
This script performs a national-level probabilistic simulation across key Italian maize-producing regions (Piemonte, Emilia Romagna, Veneto, Friuli Venezia Giulia, Lombardia) to evaluate the net economic benefits and adoption risks of the nitrogen-fixing technology. Specifically, it executes the following computational workflow:
Copula Calibration & Multivariate Modeling: Implements a calibration function (calibrate_latent_rho) to establish a latent correlation structure (targeting a correlation coefficient of 0.4) using multivariate normal distributions (MASS::mvrnorm).
Stochastic Parameter Generation: Runs Monte Carlo simulations (10,000 iterations per year) to generate correlated stochastic inputs for:
Nitrogen efficiency: modeled through a truncated normal distribution (mean = 0.30, SD = 0.05, bounds [0, 1]).
Seed price increase: modeled through a triangular distribution (minimum = 0.05, maximum = 0.16, mode = 0.10).
Economic Impact Calculation: Computes annual net economic benefits per hectare by subtracting extra seed costs (derived from regional baseline seed costs and simulated price increases) from nitrogen fertilizer cost savings.
Risk Profiles & Credible Intervals: Aggregates annual results to calculate key statistical indicators, including the probability of economic loss (Benefit < 0), the probability of failing to meet the indifference threshold (Benefit < 30 €/ha), and 95% credible intervals (2.5th and 97.5th percentiles) for each year.

Threshold_analysis.R
This script implements a deterministic partial budget evaluation framework to assess farm-level adoption dynamics and macro-sectoral economic impacts under varying seed price increase scenarios (+0%, +5%, +10%, and +16%) combined with a strict economic adoption threshold of 30 €/ha. the same framework has been run with different threshold scenario (20 €/ha and 10 €/ha)
