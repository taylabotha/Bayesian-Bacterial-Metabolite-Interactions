# Load necessary libraries for Bayesian modeling, data manipulation and visualization
library(deSolve)    # For solving differential equations (if used)
library(rstan)      # Interface to Stan for Bayesian inference
library(tidyr)      # Data tidying
library(ggplot2)    # Data visualization
library(dplyr)      # Data manipulation
library(posterior)  # Posterior analysis tools
library(stringr)    # String manipulation
library(cowplot)    # Combine multiple plots

# Set working directory for project files
setwd("/Users/tayni/Downloads/thesis proj- microbiota")

#===============================================================================
# 1. LOAD SAVED STAN FIT RESULTS
#===============================================================================
# Read in fitted Stan models for different priors and sparsity levels

# Load Uniform prior fit with 25% sparsity (path trimmed for privacy)
uniform_25 <- readRDS("......25.rds")

# Load Horseshoe prior fits, with 0%, 25%, and 50% sparsity
horseshoe_0 <- readRDS("horseshoe_fit_10x3_0%1000-500.rds")
horseshoe_25 <- readRDS("horseshoe_fit_10x3_25%1000-500.rds")
horseshoe_50 <- readRDS("horseshoe_fit_10x3_50%1000-500.rds")

# Load LASSO prior fits with 0%, 25%, and 50% sparsity
lasso_0 <- readRDS("lasso_fit_10x3_0%1000-500.rds")
lasso_25 <- readRDS("lasso_fit_10x3_25%1000-500.rds")
lasso_50 <- readRDS("lasso_fit_10x3_50%1000-500.rds")

# Load original Uniform prior fits with 0% and 50% sparsity for 10x3 community
uniform_0_correct <- readRDS("C:/Users/tayni/Downloads/thesis proj- microbiota/uniform-0%spar-10-3split-1000-500/.......rds")
uniform_50_correct <- readRDS("C:/Users/tayni/Downloads/thesis proj- microbiota/uniform-50%spar-10-3split-1000-500/......50.rds")

# Organize all fits into a single list, extracting the main fit object for each
fits <- list(
  uniform_0_correct      = uniform_0_correct[[1]],
  uniform_25             = uniform_25[[1]],
  uniform_50_correct     = uniform_50_correct[[1]],
  
  horseshoe_0            = horseshoe_0[[1]],
  horseshoe_25           = horseshoe_25[[1]],
  horseshoe_50           = horseshoe_50[[1]],
  
  lasso_0                = lasso_0[[1]],
  lasso_25               = lasso_25[[1]],
  lasso_50               = lasso_50[[1]]
)

# Create a named list for traceplot analysis with friendly labels
fits_trace <- list(
  "Uniform 0%" = fits[["uniform_0_correct"]],
  "Uniform 25%" = fits[["uniform_25"]],
  "Uniform 50%" = fits[["uniform_50_correct"]],
  
  "Horseshoe 0%" = fits[["horseshoe_0"]],
  "Horseshoe 25%" = fits[["horseshoe_25"]],
  "Horseshoe 50%" = fits[["horseshoe_50"]],
  
  "LASSO 0%" = fits[["lasso_0"]],
  "LASSO 25%" = fits[["lasso_25"]],
  "LASSO 50%" = fits[["lasso_50"]]
)

#===============================================================================
# 2. LOAD TRUE PARAMETER MATRICES FOR COVERAGE CALCULATION
#===============================================================================
truth_0  <- readRDS("truth_10x3_uni0%1000-500.rds")
truth_25 <- readRDS("truth_10x3_uni25%1000-500.rds")
truth_50 <- readRDS("truth_10x3_uni50%1000-500.rds")

# Organize true parameter matrices by sparsity level for easy access
truths <- list(
  "0" = truth_0[[1]],
  "25" = truth_25[[1]],
  "50" = truth_50[[1]]
)

#===============================================================================
# 3. CHECK FITTED PARAMETER NAMES (OPTIONAL)
#===============================================================================
# Useful to explore parameter names from the first fit
names(rstan::extract(fits[[1]], permuted = TRUE))

#===============================================================================
# 4. FUNCTION TO CALCULATE COVERAGE
#===============================================================================
# This function calculates the coverage percentage:
# the percentage of true parameter values that fall within the 95% credible intervals of the posterior.

calculate_coverage <- function(fit, true_Cij) {
  # Extract posterior samples for parameter Cij
  posterior <- rstan::extract(fit, pars = "Cij", permuted = TRUE)$Cij
  
  # Calculate 2.5% and 97.5% quantiles for each coefficient (credible interval)
  lower <- apply(posterior, c(2, 3), quantile, probs = 0.025)
  upper <- apply(posterior, c(2, 3), quantile, probs = 0.975)
  
  # Check if true values lie within the credible interval
  covered <- (true_Cij >= lower) & (true_Cij <= upper)
  
  # Return average coverage percentage
  mean(covered) * 100
}

#===============================================================================
# 5. CALCULATE COVERAGE FOR ALL MODELS
#===============================================================================
# Prepare dataframe to store coverage results for each prior and sparsity level
coverage_results <- data.frame(
  prior = rep(c("Uniform", "Horseshoe", "LASSO"), each = 3),
  sparsity = rep(c(0, 25, 50), times = 3),
  coverage = NA_real_
)

# Loop through each model, calculate coverage and store it
for (i in seq_len(nrow(coverage_results))) {
  fit <- fits[[i]]
  s <- as.character(coverage_results$sparsity[i])
  true_Cij <- truths[[s]]
  coverage_results$coverage[i] <- calculate_coverage(fit, true_Cij)
}

# Print coverage table
print(coverage_results)

#===============================================================================
# 6. PLOT COVERAGE RESULTS AS BAR CHART
#===============================================================================
coverage_plot_data <- coverage_results %>%
  mutate(
    prior = factor(as.character(prior), levels = c("Uniform", "Horseshoe", "LASSO")),
    sparsity = factor(
      as.character(sparsity),
      levels = c("0", "25", "50", "0%", "25%", "50%"),
      labels = rep(c("0%", "25%", "50%"), 2)
    )
  )

# Bar plot showing coverage by prior and sparsity
coverage_plot <- ggplot(coverage_plot_data, aes(x = sparsity, y = coverage, fill = prior)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Uniform" = "#0072B2", "Horseshoe" = "#E69F00", "LASSO" = "#009E73")) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0))) +
  labs(x = "Sparsity", y = "Coverage Percentage (%)", fill = "Prior") +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(coverage_plot)

#===============================================================================
# 7. ALTERNATIVE LINE PLOT OF COVERAGE
#===============================================================================
coverage_plot_data <- coverage_results %>%
  mutate(
    prior = factor(as.character(prior), levels = c("Uniform", "Horseshoe", "LASSO")),
    sparsity = as.numeric(gsub("%", "", as.character(sparsity)))  # Convert sparsity labels to numeric
  )

# Line plot with points showing coverage trends
coverage_plot <- ggplot(coverage_plot_data, aes(x = sparsity, y = coverage, color = prior, group = prior)) +
  geom_line(linewidth = 1, position = position_dodge(width = 1)) +
  geom_point(size = 3, position = position_dodge(width = 1)) +
  scale_color_manual(values = c("Uniform" = "#0072B2", "Horseshoe" = "#E69F00", "LASSO" = "#009E73")) +
  scale_x_continuous(breaks = c(0, 25, 50), labels = c("0%", "25%", "50%")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "Sparsity", y = "Coverage Percentage (%)", color = "Prior") +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(coverage_plot)

#===============================================================================
# 8. EFFECTIVE SAMPLE SIZE (ESS) ANALYSIS
#===============================================================================
# Extract ESS values from all fitted models for parameters Cij and sigma_ar

ESS_data <- lapply(names(fits_trace), function(model_name) {
  fit <- fits_trace[[model_name]]
  
  # Summary statistics including effective sample size
  summ <- rstan::summary(fit, pars = c("Cij", "sigma_ar"))$summary
  
  data.frame(
    model = model_name,
    parameter = rownames(summ),
    neff = summ[, "n_eff"]
  )
}) %>% bind_rows() %>%
  mutate(
    prior = case_when(
      str_detect(model, "^Uniform") ~ "Uniform",
      str_detect(model, "^Horseshoe") ~ "Horseshoe",
      str_detect(model, "^LASSO") ~ "LASSO"
    )
  )

# Plot histogram of ESS values to assess sampling efficiency
ESS_plot <- ggplot(ESS_data %>% filter(is.finite(neff)), aes(x = neff, fill = prior)) +
  geom_histogram(binwidth = 25, position = "identity", alpha = 0.5, color = "black") +
  geom_vline(xintercept = 200, linetype = "dashed", color = "black", linewidth = 0.8) + # ESS threshold
  scale_x_continuous(breaks = seq(0, 5000, 1000), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, 6000)) +
  scale_fill_manual(values = c("Uniform" = "#0072B2", "Horseshoe" = "#E69F00", "LASSO" = "#009E73")) +
  labs(x = "Effective Sample Size (ESS)", y = "Number of Parameters", fill = "Prior") +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

print(ESS_plot)

#===============================================================================
# (Additional blocks further down repeat coverage analysis and plotting for different community sizes)
#===============================================================================
