#===============================================================================
# Analysis_pipeline.R
#
# ANALYSIS PIPELINE
# Compare Uniform / Horseshoe / LASSO priors
# across 0%, 25%, 50% sparsity
#
# This script uses the files produced by the interaction script.
# It performs:
# 1. Load the 9 fitted models and shared truth/data
# 2. HMC diagnostics
# 3. Rhat / n_eff and traceplots
# 4. Posterior summaries
# 5. MAE / RMSE / correlation
# 6. Zero/non-zero recovery: sensitivity, specificity, accuracy
# 7. Credible intervals: width and coverage
# 8. Posterior predictive checks
#
# IMPORTANT:
# The same simulated data are used for all three priors within each
# sparsity level. Therefore comparisons of priors are within-sparsity
# comparisons.
#
# NOTE (what the code below actually implements): of the items listed above,
# only MAE is calculated (RMSE and correlation are not), and credible-interval
# width and coverage are not calculated; credible intervals are only shown
# graphically (see the Cij plots near the end). Section numbers in the script
# are not strictly sequential (e.g. there are two "Section 5" blocks).
#
# Inputs (all produced by the interaction/simulation script):
#   fit_<prior>_<sparsity>pct_sparsity.rds          - fitted Stan models
#   truth_Cij_<sparsity>pct_sparsity.rds            - true Cij matrices
#   data_y_obs_<sparsity>pct_sparsity.rds           - observed (simulated) data
#   posterior_Cij_<prior>_<sparsity>pct_sparsity.rds - posterior median Cij
#   y_sim_<prior>_<sparsity>pct_sparsity.rds        - posterior predictive draws
#===============================================================================

# Set working directory
# NOTE: this path is specific to the author's machine; change it before running.
setwd("/Users/tayni/Downloads/thesis proj- microbiota")

# Libraries needed
library(deSolve)    # ODE solvers (used in the simulation script; not used directly here)
library(rstan)      # Stan interface: summaries, traceplots, sampler diagnostics
library(tidyr)      # pivot_longer() for reshaping results
library(ggplot2)    # plotting
library(dplyr)      # data manipulation (bind_rows, filter, mutate, %>%)
library(posterior)  # as_draws_df() for working with posterior draws
library(stringr)    # string helpers (str_extract, str_detect)
library(bayesplot)  # (loaded but not used directly in this script)

#===============================================================================
# 0. Settings
#===============================================================================

# Your interaction script used these values
N_met <- 20                              # number of metabolites
N_bac <- 10                              # number of bacterial species
Reps  <- 6                               # number of replicates
Obs   <- 4                               # number of observed time points
delta <- 1e-6                            # small offset added to the simulated data
# (subtracted again when plotting, Section 8)
ts_sim <- seq(1e-6, 4, length.out = Obs) # observation times used in the simulation

# Use this threshold only for classifying an estimated interaction as
# practically non-zero when calculating recovery metrics.
# The true matrices contain exact zeros, so true classification is exact.
classification_threshold <- 0.1

#===============================================================================
# 1. Load the 9 fitted models and the 3 shared datasets/truths
#===============================================================================

# Named list of the nine stanfit objects: 3 priors x 3 sparsity levels.
# Names follow <prior>_<sparsity> and are parsed by the helper functions below.
fits <- list(
  uniform_0    = readRDS("fit_uniform_0pct_sparsity.rds"),
  uniform_25   = readRDS("fit_uniform_25pct_sparsity.rds"),
  uniform_50   = readRDS("fit_uniform_50pct_sparsity.rds"),
  
  horseshoe_0  = readRDS("fit_horseshoe_0pct_sparsity.rds"),
  horseshoe_25 = readRDS("fit_horseshoe_25pct_sparsity.rds"),
  horseshoe_50 = readRDS("fit_horseshoe_50pct_sparsity.rds"),
  
  lasso_0      = readRDS("fit_lasso_0pct_sparsity.rds"),
  lasso_25     = readRDS("fit_lasso_25pct_sparsity.rds"),
  lasso_50     = readRDS("fit_lasso_50pct_sparsity.rds")
)

# True Cij matrices (20 x 10: one coefficient per metabolite-bacteria pair),
# indexed by sparsity level. The saved file is a list whose first element
# is the matrix.
true_Cij <- list(
  `0`  = readRDS("truth_Cij_0pct_sparsity.rds")[[1]],
  `25` = readRDS("truth_Cij_25pct_sparsity.rds")[[1]],
  `50` = readRDS("truth_Cij_50pct_sparsity.rds")[[1]]
)

# Observed data (one simulated dataset per sparsity level, shared by all three
# priors at that sparsity level)
y_obs_all <- list(
  `0`  = readRDS("data_y_obs_0pct_sparsity.rds"),
  `25` = readRDS("data_y_obs_25pct_sparsity.rds"),
  `50` = readRDS("data_y_obs_50pct_sparsity.rds")
)

# Quick check that all nine models were loaded
cat("Number of fitted models:", length(fits), "\n")
print(names(fits))

#===============================================================================
# Helper functions
#===============================================================================

# Convert a model name such as "horseshoe_25" into prior and sparsity

# Prior name: everything before the first underscore ("horseshoe")
get_prior <- function(model_name) {
  str_extract(model_name, "^[^_]+")
}

# Sparsity level: the digits after the underscore at the end of the name (25)
get_sparsity <- function(model_name) {
  as.numeric(str_extract(model_name, "(?<=_)\\d+$"))
}

#===============================================================================
# 2. HMC diagnostics
#===============================================================================
# Hamiltonian Monte Carlo (HMC) sampler diagnostics for each of the nine models:
#   - divergent transitions (should be 0)
#   - iterations hitting the maximum tree depth (should be 0)
#   - an approximate E-BFMI (energy Bayesian fraction of missing information;
#     low values suggest the sampler explored the posterior poorly)

cat("\n============================================================\n")
cat("HMC DIAGNOSTICS\n")
cat("============================================================\n")

hmc_diagnostics <- lapply(names(fits), function(model_name) {
  
  cat("\n==== Diagnostics for", model_name, "====\n")
  
  fit <- fits[[model_name]]
  
  # Print Stan's diagnostic messages
  rstan::check_hmc_diagnostics(fit)
  
  # Also save the key diagnostic counts
  # (one matrix of sampler statistics per chain, post-warm-up only)
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  
  # Total number of divergent transitions across all chains
  divergences <- sum(
    sapply(sampler_params, function(x) {
      if ("divergent__" %in% colnames(x)) sum(x[, "divergent__"])
      else 0
    })
  )
  
  # Total number of iterations that reached the maximum tree depth.
  # NOTE: assumes Stan's default max_treedepth of 10.
  max_treedepth <- sum(
    sapply(sampler_params, function(x) {
      if ("treedepth__" %in% colnames(x)) {
        max_td <- max(x[, "treedepth__"])
        if (max_td >= 10) sum(x[, "treedepth__"] >= 10) else 0
      } else {
        0
      }
    })
  )
  
  # Energy values from every chain, combined into one vector
  energy <- lapply(sampler_params, function(x) {
    if ("energy__" %in% colnames(x)) x[, "energy__"] else NULL
  })
  
  energy <- unlist(energy)
  
  # Approximate E-BFMI: mean squared successive change in energy divided by
  # the variance of the energy. This is an approximation because the chains are
  # concatenated, so the differences at the boundary between chains are
  # included. Very low values (commonly below about 0.2-0.3) are a warning sign.
  bfmi <- if (length(energy) > 1) {
    mean(diff(energy)^2) / var(energy)
  } else {
    NA_real_
  }
  
  data.frame(
    model = model_name,
    divergences = divergences,
    max_treedepth_hits = max_treedepth,
    E_BFMI_approx = bfmi
  )
})

# Combine into one table with one row per model
hmc_diagnostics_df <- bind_rows(hmc_diagnostics)

print(hmc_diagnostics_df)

#===============================================================================
# 3. Traceplots
#===============================================================================

# Cij has 200 parameters, so plotting every Cij trace at once is very crowded.
# Here we plot sigma_ar and a few representative Cij parameters.
# NOTE: the code below plots only the four Cij parameters listed (sigma_ar is
# not included in `pars`).

# One traceplot (all chains overlaid) per model
traceplots <- lapply(names(fits), function(model_name) {
  
  fit <- fits[[model_name]]
  
  traceplot(
    fit,
    pars = c(
      "Cij[1,1]",
      "Cij[5,6]",
      "Cij[15,4]",
      "Cij[18,3]"
    )
  ) +
    ggtitle(paste("Traceplot:", model_name))
})

# Name the list elements by model so individual plots can be retrieved by name
names(traceplots) <- names(fits)

# Print examples: one model per prior
print(traceplots[["uniform_0"]])
print(traceplots[["horseshoe_25"]])
print(traceplots[["lasso_50"]])

#===============================================================================
# 4. Rhat, n_eff and posterior summaries
#===============================================================================

# Posterior summaries for sigma_ar and all Cij parameters
# (mean, sd, quantiles, n_eff and Rhat for each parameter)
params_list <- lapply(fits, function(f) {
  
  rstan::summary(
    f,
    pars = c("sigma_ar", "Cij")
  )$summary
  
})

names(params_list) <- names(fits)

# Example
print(params_list[["uniform_0"]])

# Extract Rhat and n_eff for all parameters in all models.
# bind_rows(.id = "model") records the model name in a "model" column.
diagnostics_df <- lapply(fits, function(f) {
  
  summ <- rstan::summary(f)$summary
  
  data.frame(
    parameter = rownames(summ),
    Rhat = summ[, "Rhat"],
    neff = summ[, "n_eff"]
  )
  
}) %>%
  bind_rows(.id = "model")

# Keep only the interaction parameters and sigma_ar.
# NOTE: this filter assumes the Stan model declares parameters named
# Cij_pos[...] and Cij_neg[...] (the sampled coefficients) alongside the
# derived Cij matrix. If those names do not exist, only sigma_ar is kept.
selected_diagnostics <- diagnostics_df %>%
  filter(
    str_detect(parameter, "Cij_pos\\[") |
      str_detect(parameter, "Cij_neg\\[") |
      parameter == "sigma_ar"
  )

# Print parameters with Rhat above 1.01 (a common convergence threshold)
poor_rhat <- selected_diagnostics %>%
  filter(is.finite(Rhat), Rhat > 1.01)

cat("\nParameters with Rhat > 1.01:\n")
print(poor_rhat)

# Rhat histogram, one panel per model; dashed line marks the 1.01 threshold
Rhat_plot <- ggplot(
  selected_diagnostics %>% filter(is.finite(Rhat)),
  aes(x = Rhat)
) +
  geom_histogram(binwidth = 0.0005) +
  geom_vline(xintercept = 1.01, linetype = "dashed") +
  facet_wrap(~ model, scales = "free_x", ncol = 3) +
  labs(
    x = "Rhat",
    y = "Number of parameters",
    title = "Rhat diagnostics across all fitted models"
  ) +
  theme_minimal()

print(Rhat_plot)

# n_eff histogram, one panel per model.
# NOTE: scale_y_continuous(limits = c(1, 10)) removes any histogram bar taller
# than 10 parameters (and gives a warning), which can hide bins with many
# parameters. To zoom in without dropping data, use
# coord_cartesian(ylim = c(0, 10)) instead.
N_eff_plot <- ggplot(
  selected_diagnostics %>% filter(is.finite(neff)),
  aes(x = neff)
) +
  geom_histogram(binwidth = 25) +
  facet_wrap(~ model, scales = "free_x") + scale_y_continuous(breaks = 1:10, limits = c(1,10))+
  labs(
    x = "Effective sample size",
    y = "Number of parameters",
    title = "Effective sample size across all fitted models"
  ) +
  theme_minimal()

print(N_eff_plot)

#-------------------------------------------------------------------------------
# Sanity check: are the Uniform 0% and Uniform 50% fits different?
#-------------------------------------------------------------------------------
# These models were fitted to different data, so their posteriors should differ.
# If the two comparisons below report equality, the same fit file was probably
# saved under two names.

# Extract posterior medians directly from each Stan fit
# ("50%" is the median column of rstan's summary matrix)
uniform_0_medians <- rstan::summary(
  fits[["uniform_0"]],
  pars = "Cij"
)$summary[, "50%"]

uniform_50_medians <- rstan::summary(
  fits[["uniform_50"]],
  pars = "Cij"
)$summary[, "50%"]

# Compare them (TRUE would mean the medians are identical up to tolerance)
all.equal(
  uniform_0_medians,
  uniform_50_medians
)

# Compare the full sets of posterior draws as well
draws_0 <- rstan::extract(
  fits[["uniform_0"]],
  pars = "Cij"
)$Cij

draws_50 <- rstan::extract(
  fits[["uniform_50"]],
  pars = "Cij"
)$Cij

identical(draws_0, draws_50)

#===============================================================================
# 5. Posterior estimates versus true Cij 50% sparsity
#===============================================================================
# NOTE: despite the heading, this block loads the posterior medians for all
# nine models and calculates MAE for every model, not only the 50% sparsity ones.

# We use the saved posterior median Cij values produced by the interaction
# script. These are directly comparable with the true Cij matrix.

posterior_Cij <- lapply(names(fits), function(model_name) {
  
  # (this value is not used inside this function)
  sparsity <- get_sparsity(model_name)
  
  readRDS(
    paste0(
      "posterior_Cij_",
      model_name,
      "pct_sparsity.rds"
    )
  )
  
})

names(posterior_Cij) <- names(fits)

# IMPORTANT:
# The filenames are "posterior_Cij_uniform_0pct_sparsity.rds", etc.
# The model_name already contains the prior and sparsity number.

# Calculate MAE, RMSE and correlation for every model
# NOTE: only MAE is calculated below. Mean absolute error = average absolute
# difference between the true and estimated Cij over all 200 coefficients.
accuracy_results <- lapply(names(fits), function(model_name) {
  
  sparsity <- get_sparsity(model_name)
  
  # Flatten the true and estimated matrices into vectors for comparison
  true_vals <- as.vector(true_Cij[[as.character(sparsity)]])
  est_vals <- as.vector(posterior_Cij[[model_name]])
  
  data.frame(
    model = model_name,
    prior = get_prior(model_name),
    sparsity = sparsity,
    MAE = mean(abs(true_vals - est_vals))
  )
  
}) %>%
  bind_rows()

print(accuracy_results)




#===============================================================================
# 5b. TRUE Cij VS POSTERIOR MEDIAN
#     20 x 10 community, 50% sparsity
#     Compare Uniform, Horseshoe and LASSO
#===============================================================================

# True Cij values for the 50% sparse community
true_vals_50 <- as.vector(true_Cij[["50"]])


#-------------------------------------------------------------------------------
# Create dataframe for the three priors
#-------------------------------------------------------------------------------
# Long data frame with one row per coefficient per prior: the true value
# (the same for all priors) next to each prior's posterior median.

Cij_recovery_50 <- bind_rows(
  
  data.frame(
    prior = "Uniform",
    true_Cij = true_vals_50,
    posterior_median = as.vector(
      posterior_Cij[["uniform_50"]]
    )
  ),
  
  data.frame(
    prior = "Horseshoe",
    true_Cij = true_vals_50,
    posterior_median = as.vector(
      posterior_Cij[["horseshoe_50"]]
    )
  ),
  
  data.frame(
    prior = "LASSO",
    true_Cij = true_vals_50,
    posterior_median = as.vector(
      posterior_Cij[["lasso_50"]]
    )
  )
)


#-------------------------------------------------------------------------------
# Identify whether each TRUE interaction is zero or non-zero
#-------------------------------------------------------------------------------

Cij_recovery_50 <- Cij_recovery_50 %>%
  mutate(
    interaction_type = ifelse(
      true_Cij == 0,
      "True zero",
      "True non-zero"
    )
  )


# Make sure panels appear in desired order
Cij_recovery_50$prior <- factor(
  Cij_recovery_50$prior,
  levels = c(
    "Uniform",
    "Horseshoe",
    "LASSO"
  )
)


#-------------------------------------------------------------------------------
# Check the dataframe
#-------------------------------------------------------------------------------

print(head(Cij_recovery_50))

# Number of true zero / true non-zero coefficients for each prior
table(
  Cij_recovery_50$prior,
  Cij_recovery_50$interaction_type
)


#===============================================================================
# PLOT
#===============================================================================
# Scatter plot of true Cij (x) against posterior median (y), one panel per
# prior. Points on the dashed 1:1 line are perfectly recovered.

Cij_recovery_plot <- ggplot(
  Cij_recovery_50,
  aes(
    x = true_Cij,
    y = posterior_median,
    shape = interaction_type,
    colour = interaction_type
  )
) +
  
  geom_point(
    size = 2,
    alpha = 0.7
  ) +
  
  # Perfect recovery line
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.8,
    colour = "black"
  ) +
  
  # Three priors
  facet_wrap(
    ~ prior,
    nrow = 1
  ) +
  
  # Same axes for all three panels
  coord_equal(
    xlim = c(-1, 1),
    ylim = c(-1, 1)
  ) +
  
  scale_x_continuous(
    breaks = seq(-1, 1, by = 0.5)
  ) +
  
  scale_y_continuous(
    breaks = seq(-1, 1, by = 0.5)
  ) +
  
  # Colour-blind-safe colours
  scale_colour_manual(
    values = c(
      "True non-zero" = "#0072B2",
      "True zero" = "#E69F00"
    )
  ) +
  
  # Also distinguish by shape
  scale_shape_manual(
    values = c(
      "True non-zero" = 16,
      "True zero" = 17
    )
  ) +
  
  labs(
    x = expression("True interaction coefficient (" * C[ij] * ")"),
    y = expression("Posterior median (" * C[ij] * ")"),
    colour = "True interaction",
    shape = "True interaction"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title.x = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.position = "bottom",
    
    panel.grid.minor = element_blank()
  )

print(Cij_recovery_plot)













#===============================================================================
# 6. Zero / non-zero interaction recovery
#===============================================================================

# True values:
#   zero       = exactly 0
#   non-zero   = anything else
#
# Estimated values:
#   estimated non-zero if |posterior median| > classification_threshold
#
# This allows sensitivity, specificity and accuracy to be calculated.
#
# Metric definitions (TP/TN/FP/FN = true/false positives/negatives, where
# "positive" means "non-zero interaction"):
#   sensitivity = TP / (TP + FN)   proportion of real interactions detected
#   specificity = TN / (TN + FP)   proportion of true zeros identified as zero
#   accuracy    = (TP + TN) / total number of coefficients
#   precision   = TP / (TP + FP)   proportion of detected interactions that are real
# A metric is NA when its denominator is zero (e.g. specificity at 0% sparsity,
# where there are no true zeros).

recovery_results <- lapply(names(fits), function(model_name) {
  
  sparsity <- get_sparsity(model_name)
  
  true_vals <- as.vector(true_Cij[[as.character(sparsity)]])
  est_vals <- as.vector(posterior_Cij[[model_name]])
  
  # Classify each coefficient as non-zero or not
  true_nonzero <- abs(true_vals) > 0
  estimated_nonzero <- abs(est_vals) > classification_threshold
  
  # Confusion matrix counts
  TP <- sum(true_nonzero & estimated_nonzero)     # real interaction, detected
  TN <- sum(!true_nonzero & !estimated_nonzero)   # true zero, estimated as zero
  FP <- sum(!true_nonzero & estimated_nonzero)    # true zero, wrongly detected
  FN <- sum(true_nonzero & !estimated_nonzero)    # real interaction, missed
  
  sensitivity <- ifelse(
    (TP + FN) > 0,
    TP / (TP + FN),
    NA
  )
  
  specificity <- ifelse(
    (TN + FP) > 0,
    TN / (TN + FP),
    NA
  )
  
  accuracy <- (TP + TN) / length(true_vals)
  
  precision <- ifelse(
    (TP + FP) > 0,
    TP / (TP + FP),
    NA
  )
  
  data.frame(
    model = model_name,
    prior = get_prior(model_name),
    sparsity = sparsity,
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    sensitivity = sensitivity,
    specificity = specificity,
    accuracy = accuracy,
    precision = precision
  )
  
}) %>%
  bind_rows()

print(recovery_results)





#===============================================================================
# FIGURE 3: MAE VS SPARSITY
#===============================================================================
# Line plot of MAE against sparsity level, one line per prior.
# (Prior names are lower case in the data, so the colour/label vectors are
# keyed by the lower-case names and relabelled for display.)

MAE_plot <- ggplot(
  accuracy_results,
  aes(
    x = sparsity,
    y = MAE,
    colour = prior,
    group = prior
  )
) +
  
  geom_line(linewidth = 0.9) +
  
  geom_point(size = 3) +
  
  scale_x_continuous(
    breaks = c(0, 25, 50),
    labels = c("0%", "25%", "50%")
  ) +
  
  scale_colour_manual(
    values = c(
      "uniform" = "#0072B2",
      "horseshoe" = "#D55E00",
      "lasso" = "#009E73"
    ),
    labels = c(
      "uniform" = "Uniform",
      "horseshoe" = "Horseshoe",
      "lasso" = "LASSO"
    )
  ) +
  
  labs(
    x = "Sparsity",
    y = "Mean absolute error (MAE)",
    colour = "Prior"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(MAE_plot)


#===============================================================================
# FIGURE 4: SENSITIVITY AND SPECIFICITY
#===============================================================================

# Convert sensitivity and specificity into long format
# (one row per model-metric pair, so ggplot can facet by metric)
sensitivity_specificity_df <- recovery_results %>%
  select(
    prior,
    sparsity,
    sensitivity,
    specificity
  ) %>%
  
  pivot_longer(
    cols = c(sensitivity, specificity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  
  mutate(
    metric = factor(
      metric,
      levels = c("sensitivity", "specificity"),
      labels = c("Sensitivity", "Specificity")
    )
  )


#-------------------------------------------------------------------------------
# Create two-panel graph
#-------------------------------------------------------------------------------
# One line per prior across the three sparsity levels, with separate panels
# for sensitivity and specificity. NA values (e.g. specificity at 0% sparsity)
# are dropped from the lines and points by na.rm = TRUE.

Sensitivity_Specificity_plot <- ggplot(
  sensitivity_specificity_df,
  aes(
    x = sparsity,
    y = value,
    colour = prior,
    group = prior
  )
) +
  
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  
  geom_point(size = 3, na.rm = TRUE) +
  
  facet_wrap(
    ~ metric,
    ncol = 2,
    scales = "free_x"
  )+
  
  scale_x_continuous(
    breaks = c(0, 25, 50),
    labels = c("0%", "25%", "50%")
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  
  scale_colour_manual(
    values = c(
      "uniform" = "#0072B2",
      "horseshoe" = "#D55E00",
      "lasso" = "#009E73"
    ),
    labels = c(
      "uniform" = "Uniform",
      "horseshoe" = "Horseshoe",
      "lasso" = "LASSO"
    )
  ) +
  
  labs(
    x = "Sparsity",
    y = "Proportion correctly classified",
    colour = "Prior"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(Sensitivity_Specificity_plot)





#===============================================================================
# TRUE-ZERO VS TRUE-NON-ZERO MAE
#===============================================================================
# Split the error into two parts:
#   MAE_zero    - average |estimate| over coefficients whose true value is 0
#                 (any non-zero estimate is an error)
#   MAE_nonzero - average |true - estimate| over the truly non-zero coefficients
# NA is returned if a group is empty (e.g. no true zeros at 0% sparsity).

recovery_MAE <- lapply(names(posterior_Cij), function(model_name) {
  
  sparsity <- get_sparsity(model_name)
  
  true_vals <- as.vector(
    true_Cij[[as.character(sparsity)]]
  )
  
  est_vals <- as.vector(
    posterior_Cij[[model_name]]
  )
  
  # Logical indices of the true zeros and true non-zeros
  zero_idx <- true_vals == 0
  nonzero_idx <- true_vals != 0
  
  data.frame(
    model = model_name,
    prior = get_prior(model_name),
    sparsity = sparsity,
    
    MAE_zero = if (any(zero_idx)) {
      mean(abs(est_vals[zero_idx]))
    } else {
      NA_real_
    },
    
    MAE_nonzero = if (any(nonzero_idx)) {
      mean(abs(
        true_vals[nonzero_idx] -
          est_vals[nonzero_idx]
      ))
    } else {
      NA_real_
    }
  )
  
}) %>%
  bind_rows()


# Print results
print(recovery_MAE)


#===============================================================================
# FIGURE 5: ZERO VS NON-ZERO COEFFICIENT RECOVERY
#===============================================================================

# Convert MAE results to long format
MAE_comparison_df <- recovery_MAE %>%
  
  select(
    prior,
    sparsity,
    MAE_zero,
    MAE_nonzero
  ) %>%
  
  pivot_longer(
    cols = c(MAE_zero, MAE_nonzero),
    names_to = "metric",
    values_to = "MAE"
  ) %>%
  
  mutate(
    metric = factor(
      metric,
      levels = c("MAE_zero", "MAE_nonzero"),
      labels = c(
        "True-zero coefficients",
        "True-non-zero coefficients"
      )
    )
  )


#-------------------------------------------------------------------------------
# Create two-panel graph
#-------------------------------------------------------------------------------
# MAE against sparsity for true-zero (left) and true-non-zero (right)
# coefficients, one line per prior.

MAE_recovery_plot <- ggplot(
  MAE_comparison_df,
  aes(
    x = sparsity,
    y = MAE,
    colour = prior,
    group = prior
  )
) +
  
  geom_line(
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  
  geom_point(
    size = 3,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ metric,
    ncol = 2
  ) +
  
  scale_x_continuous(
    breaks = c(0, 25, 50),
    labels = c("0%", "25%", "50%")
  ) +
  
  scale_colour_manual(
    values = c(
      "uniform" = "#0072B2",
      "horseshoe" = "#D55E00",
      "lasso" = "#009E73"
    ),
    labels = c(
      "uniform" = "Uniform",
      "horseshoe" = "Horseshoe",
      "lasso" = "LASSO"
    )
  ) +
  
  labs(
    x = "Sparsity",
    y = "Mean absolute error (MAE)",
    colour = "Prior"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(MAE_recovery_plot)



#===============================================================================
# 6c. Accuracy plot
#===============================================================================

# Accuracy is the proportion of Cij interactions that were correctly
# classified as zero or non-zero.
#
# No error bars are shown here because there is only one accuracy estimate
# for each prior × sparsity combination. The credible intervals for the
# individual Cij coefficients are shown separately in the Cij plots.
#
# NOTE: this figure uses a different colour for each prior than the other
# figures (e.g. Horseshoe is blue here but orange elsewhere). Consider using
# one palette throughout the thesis.


Accuracy_plot <- ggplot(
  recovery_results,
  aes(
    x = prior,
    y = accuracy,
    fill = prior
  )
) +
  
  # Accuracy bars
  geom_col(
    width = 0.7,
    colour = "black",
    linewidth = 0.25
  ) +
  
  # Colour-blind-friendly Okabe-Ito palette
  scale_fill_manual(
    values = c(
      "horseshoe" = "#0072B2",  # blue
      "lasso"     = "#E69F00",  # orange
      "uniform"   = "#009E73"   # green
    ),
    labels = c(
      "horseshoe" = "Horseshoe",
      "lasso"     = "LASSO",
      "uniform"   = "Uniform"
    )
  ) +
  
  # One panel for each sparsity level
  facet_wrap(
    ~ sparsity,
    nrow = 1,
    labeller = labeller(
      sparsity = function(x) {
        paste0("Sparsity = ", x, "%")
      }
    )
  ) +
  
  # Accuracy ranges from 0 to 1
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    x = "Prior",
    y = "Classification accuracy"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    
    # Make sparsity headings prominent
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    
    # Prior names underneath bars
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    # No minor gridlines
    panel.grid.minor = element_blank(),
    
    # Slightly more space between panels
    panel.spacing = unit(1, "lines"),
    
    # No legend because the prior is already labelled on the x-axis
    legend.position = "none"
  )

print(Accuracy_plot)


#===============================================================================
# 6d. Number of detected non-zero interactions
#===============================================================================
# For each model: the number of truly non-zero coefficients compared with the
# number the model estimates as non-zero (|posterior median| > threshold).

detected_nonzero <- lapply(names(fits), function(model_name) {
  
  sparsity <- get_sparsity(model_name)
  
  true_vals <- as.vector(true_Cij[[as.character(sparsity)]])
  est_vals <- as.vector(posterior_Cij[[model_name]])
  
  data.frame(
    model = model_name,
    prior = get_prior(model_name),
    sparsity = sparsity,
    true_nonzero = sum(abs(true_vals) > 0),
    estimated_nonzero = sum(abs(est_vals) > classification_threshold)
  )
  
}) %>%
  bind_rows()

print(detected_nonzero)


#===============================================================================
# 7c. Posterior distributions for selected Cij parameters
#===============================================================================

# Rather than plotting all 200 parameters, choose representative parameters:
# one true zero/non-zero consumed coefficient and one produced coefficient.
#
# Change these indices if you want to inspect other interactions.

selected_parameters <- c(
  "Cij[1,1]",
  "Cij[5,6]",
  "Cij[15,4]",
  "Cij[20,3]"
)

# Collect the posterior draws of the selected coefficients from every model
posterior_density_data <- lapply(names(fits), function(model_name) {
  
  fit <- fits[[model_name]]
  sparsity <- get_sparsity(model_name)
  
  # All posterior draws as a data frame (one column per parameter)
  draws <- as_draws_df(fit)
  
  # Only keep selected parameters that exist in this fit
  available_parameters <- intersect(
    selected_parameters,
    names(draws)
  )
  
  if (length(available_parameters) == 0) {
    return(NULL)
  }
  
  # Reshape to long format: one row per draw per parameter
  draws %>%
    select(all_of(available_parameters)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      model = model_name,
      prior = get_prior(model_name),
      sparsity = sparsity
    )
  
}) %>%
  bind_rows()

# Add true values
# (row and column indices are parsed from the parameter name, e.g.
# "Cij[15,4]" -> row 15, column 4, then looked up in the true matrix)
posterior_density_data <- posterior_density_data %>%
  mutate(
    row_idx = as.integer(str_extract(parameter, "(?<=\\[)\\d+")),
    col_idx = as.integer(str_extract(parameter, "(?<=,)\\d+"))
  ) %>%
  rowwise() %>%
  mutate(
    true_value = true_Cij[[as.character(sparsity)]][row_idx, col_idx]
  ) %>%
  ungroup()

# Posterior density grid: rows = parameters, columns = models.
# Dashed vertical line = true value of the coefficient.
Density_plot <- ggplot(
  posterior_density_data,
  aes(x = value)
) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(
    aes(xintercept = true_value),
    linetype = "dashed",
    linewidth = 1
  ) +
  facet_grid(
    parameter ~ model,
    scales = "free"
  ) +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distributions for selected Cij parameters"
  ) +
  theme_minimal()

print(Density_plot)

#===============================================================================
# 8. Posterior predictive checks
#===============================================================================

# Function to calculate and plot the posterior predictive mean and 95%
# posterior predictive interval against the observed data.
#
# Argument: model_name, e.g. "uniform_25"
# Returns:  a list with two ggplot objects (metabolites, bacteria); both plots
#           are also printed.

make_ppc_plot <- function(model_name) {
  
  sparsity <- get_sparsity(model_name)
  
  # Load posterior predictive samples saved by interaction script
  y_sim <- readRDS(
    paste0(
      "y_sim_",
      model_name,
      "pct_sparsity.rds"
    )
  )
  
  # Observed data
  # (the first index of the saved array is dropped, leaving
  # replicate x time x variable)
  y_obs <- y_obs_all[[as.character(sparsity)]][1,,,]
  
  # Remove the delta added during simulation for plotting
  y_obs_corrected <- y_obs - delta
  
  # y_sim dimensions:
  # posterior draw x replicate x time x variable
  
  # Mean over replicates for each posterior draw
  y_sim_mean_reps <- apply(
    y_sim,
    c(1, 3, 4),
    mean
  )
  
  # Mean over posterior draws
  y_sim_mean <- apply(
    y_sim_mean_reps,
    c(2, 3),
    mean
  )
  
  # 95% posterior predictive interval
  y_sim_lower <- apply(
    y_sim_mean_reps,
    c(2, 3),
    quantile,
    probs = 0.025
  )
  
  y_sim_upper <- apply(
    y_sim_mean_reps,
    c(2, 3),
    quantile,
    probs = 0.975
  )
  
  # Observed mean and SD across replicates
  y_obs_mean <- apply(
    y_obs_corrected,
    c(2, 3),
    mean
  )
  
  y_obs_sd <- apply(
    y_obs_corrected,
    c(2, 3),
    sd
  )
  
  # Number of time points and number of variables (metabolites + bacteria).
  # NOTE: naming a variable `T` masks R's shorthand for TRUE; consider
  # renaming it (e.g. n_time).
  T <- dim(y_sim)[3]
  L <- dim(y_sim)[4]
  
  # Data frame with one row per time point per variable. expand.grid varies
  # time fastest, which matches the column-major order of as.vector() applied
  # to the time x variable matrices below.
  df_across <- expand.grid(
    time = ts_sim[1:T],
    variable = seq_len(L)
  )
  
  df_across$obs <- as.vector(y_obs_mean)       # observed mean over replicates
  df_across$obs_sd <- as.vector(y_obs_sd)      # observed SD over replicates
  df_across$mean <- as.vector(y_sim_mean)      # posterior predictive mean
  df_across$lower <- as.vector(y_sim_lower)    # 2.5% posterior predictive bound
  df_across$upper <- as.vector(y_sim_upper)    # 97.5% posterior predictive bound
  
  df_across$variable_name <- paste0(
    "Variable ",
    df_across$variable
  )
  
  # Metabolites (the first N_met variables)
  df_metabolites <- subset(
    df_across,
    variable <= N_met
  )
  
  # Bacteria (the remaining variables)
  df_bacteria <- subset(
    df_across,
    variable > N_met
  )
  
  # Metabolite plot
  # Ribbon = 95% posterior predictive interval, solid line = posterior
  # predictive mean, dashed line = observed mean, error bars = observed
  # mean +/- 1 SD across replicates.
  p_met <- ggplot(
    df_metabolites,
    aes(x = time)
  ) +
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper
      ),
      alpha = 0.3
    ) +
    geom_line(
      aes(y = mean),
      linewidth = 1
    ) +
    geom_line(
      aes(y = obs),
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_errorbar(
      aes(
        ymin = obs - obs_sd,
        ymax = obs + obs_sd
      ),
      width = 0.2
    ) +
    facet_wrap(
      ~ variable_name,
      scales = "free_y",
      ncol = 5
    ) +
    labs(
      x = "Time",
      y = "Metabolite abundance",
      title = paste(
        "Posterior predictive check:",
        model_name
      ),
      subtitle = "Solid = posterior mean, dashed = observed mean, ribbon = 95% posterior interval"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 10)
    )
  
  print(p_met)
  
  # Bacteria plot (same layout and elements as the metabolite plot)
  p_bac <- ggplot(
    df_bacteria,
    aes(x = time)
  ) +
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper
      ),
      alpha = 0.3
    ) +
    geom_line(
      aes(y = mean),
      linewidth = 1
    ) +
    geom_line(
      aes(y = obs),
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_errorbar(
      aes(
        ymin = obs - obs_sd,
        ymax = obs + obs_sd
      ),
      width = 0.2
    ) +
    facet_wrap(
      ~ variable_name,
      scales = "free_y",
      ncol = 5
    ) +
    labs(
      x = "Time",
      y = "Bacteria abundance",
      title = paste(
        "Posterior predictive check:",
        model_name
      ),
      subtitle = "Solid = posterior mean, dashed = observed mean, ribbon = 95% posterior interval"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 10)
    )
  
  print(p_bac)
  
  return(
    list(
      metabolites = p_met,
      bacteria = p_bac
    )
  )
}

# IMPORTANT:
# There are 9 models, so do not print all 9 huge PPC figures immediately.
# Start with one model from each prior at 25% sparsity.
#
# If these look correct, run the other models below.

ppc_uniform_25 <- make_ppc_plot("uniform_25")
ppc_horseshoe_25 <- make_ppc_plot("horseshoe_25")
ppc_lasso_25 <- make_ppc_plot("lasso_25")

# To make PPC plots for every model:
# ppc_all <- lapply(names(fits), make_ppc_plot)


#===============================================================================
# Posterior mean and 95% credible intervals for Cij
#===============================================================================
# The following four blocks make the same plot for four selected models
# (Uniform 25%, Horseshoe 50%, LASSO 25%, LASSO 50%). For each model, 50
# randomly chosen Cij coefficients are shown with their point estimate (black
# dot) and 95% credible interval (red bar); the black vertical line at 0 marks
# "no interaction". An interval that excludes 0 indicates a detected interaction.
#
# NOTE: the same code is repeated four times, changing only the model name.
# It could be wrapped in a function taking the model name.

#-------------------------------------------------------------------------------
# Uniform prior, 25% sparsity
#-------------------------------------------------------------------------------

fit <- fits[["uniform_25"]]

# Summary matrix for all Cij parameters (mean, quantiles, ...)
Cij_summary <- rstan::summary(
  fit,
  pars = "Cij"
)$summary

# NOTE: rstan's summary matrix has no column called "median"; the median is
# stored in the "50%" column. As written, `Cij_summary[, "median"]` is likely to
# fail with "subscript out of bounds". Use Cij_summary[, "50%"] instead.
# NOTE: this block uses the posterior MEDIAN, while the three blocks below use
# the posterior MEAN. The plot legend comments still say "Posterior mean".
Cij_df <- data.frame(
  parameter = rownames(Cij_summary),
  median = Cij_summary[, "median"],  #changed from mean to median???
  lower = Cij_summary[, "2.5%"],
  upper = Cij_summary[, "97.5%"]
)

# Randomly select 50 interactions
# (fixed seed so the same 50 are selected every time)
set.seed(42)

Cij_df <- Cij_df %>%
  slice_sample(n = 50)

# Give the selected interactions numbers 1–50 for the y-axis
Cij_df$parameter_number <- seq_len(nrow(Cij_df))

Cij_plot <- ggplot(
  Cij_df,
  aes(x = median, y = parameter_number)  #changed x from mean to median
) +
  # 95% credible interval = RED
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    height = 0.2,
    color = "red"
  ) +
  # Posterior mean = BLACK DOT
  geom_point(
    size = 2,
    color = "black"
  ) +
  # Zero interaction = BLACK vertical line
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    linewidth = 1,
    color = "black"
  ) +
  # Label each row of the plot with its parameter name (e.g. "Cij[3,7]")
  scale_y_continuous(
    breaks = Cij_df$parameter_number,
    labels = Cij_df$parameter
  ) +
  labs(
    x = "Interaction coefficient (Cij)",
    y = "Metabolite-bacteria pair"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 7)
  )

print(Cij_plot)

#-------------------------------------------------------------------------------
# Horseshoe prior, 50% sparsity
#-------------------------------------------------------------------------------
fit <- fits[["horseshoe_50"]]

Cij_summary <- rstan::summary(
  fit,
  pars = "Cij"
)$summary

# Uses the posterior mean (not the median)
Cij_df <- data.frame(
  parameter = rownames(Cij_summary),
  mean = Cij_summary[, "mean"],
  lower = Cij_summary[, "2.5%"],
  upper = Cij_summary[, "97.5%"]
)

# Randomly select 50 interactions
set.seed(42)

Cij_df <- Cij_df %>%
  slice_sample(n = 50)

# Give the selected interactions numbers 1–50 for the y-axis
Cij_df$parameter_number <- seq_len(nrow(Cij_df))

Cij_plot <- ggplot(
  Cij_df,
  aes(x = mean, y = parameter_number)
) +
  # 95% credible interval = RED
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    height = 0.2,
    color = "red"
  ) +
  # Posterior mean = BLACK DOT
  geom_point(
    size = 2,
    color = "black"
  ) +
  # Zero interaction = BLACK vertical line
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    linewidth = 1,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = Cij_df$parameter_number,
    labels = Cij_df$parameter
  ) +
  labs(
    x = "Interaction coefficient (Cij)",
    y = "Metabolite-bacteria pair"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 7)
  )

print(Cij_plot)

#-------------------------------------------------------------------------------
# LASSO prior, 25% sparsity
#-------------------------------------------------------------------------------
fit <- fits[["lasso_25"]]

Cij_summary <- rstan::summary(
  fit,
  pars = "Cij"
)$summary

Cij_df <- data.frame(
  parameter = rownames(Cij_summary),
  mean = Cij_summary[, "mean"],
  lower = Cij_summary[, "2.5%"],
  upper = Cij_summary[, "97.5%"]
)

# Randomly select 50 interactions
set.seed(42)

Cij_df <- Cij_df %>%
  slice_sample(n = 50)

# Give the selected interactions numbers 1–50 for the y-axis
Cij_df$parameter_number <- seq_len(nrow(Cij_df))

Cij_plot <- ggplot(
  Cij_df,
  aes(x = mean, y = parameter_number)
) +
  # 95% credible interval = RED
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    height = 0.2,
    color = "red"
  ) +
  # Posterior mean = BLACK DOT
  geom_point(
    size = 2,
    color = "black"
  ) +
  # Zero interaction = BLACK vertical line
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    linewidth = 1,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = Cij_df$parameter_number,
    labels = Cij_df$parameter
  ) +
  labs(
    x = "Interaction coefficient (Cij)",
    y = "Metabolite-bacteria pair"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 7)
  )

print(Cij_plot)

#-------------------------------------------------------------------------------
# LASSO prior, 50% sparsity
#-------------------------------------------------------------------------------
fit <- fits[["lasso_50"]]

Cij_summary <- rstan::summary(
  fit,
  pars = "Cij"
)$summary

Cij_df <- data.frame(
  parameter = rownames(Cij_summary),
  mean = Cij_summary[, "mean"],
  lower = Cij_summary[, "2.5%"],
  upper = Cij_summary[, "97.5%"]
)

# Randomly select 50 interactions
set.seed(42)

Cij_df <- Cij_df %>%
  slice_sample(n = 50)

# Give the selected interactions numbers 1–50 for the y-axis
Cij_df$parameter_number <- seq_len(nrow(Cij_df))

Cij_plot <- ggplot(
  Cij_df,
  aes(x = mean, y = parameter_number)
) +
  # 95% credible interval = RED
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    height = 0.2,
    color = "red"
  ) +
  # Posterior mean = BLACK DOT
  geom_point(
    size = 2,
    color = "black"
  ) +
  # Zero interaction = BLACK vertical line
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    linewidth = 1,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = Cij_df$parameter_number,
    labels = Cij_df$parameter
  ) +
  labs(
    x = "Interaction coefficient (Cij)",
    y = "Metabolite-bacteria pair"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 7)
  )

print(Cij_plot)


#===============================================================================
# 9. Final combined results table
#===============================================================================
# Combine the MAE (accuracy_results) and the zero/non-zero classification
# accuracy (recovery_results) into one table, sorted by sparsity then prior.

final_results <- accuracy_results %>%
  left_join(
    recovery_results %>%
      select(
        model,
        accuracy
      ),
    by = "model"
  ) %>%
  arrange(sparsity, prior)

print(final_results)

# Compact table of classification accuracy only
accuracy_table <- recovery_results %>%
  select(prior, sparsity, accuracy) %>%
  arrange(sparsity, prior)

print(accuracy_table)


#===============================================================================
# 5. Posterior estimates versus true Cij all sparsity and prior levels
#zero vs non-zero graphs
#===============================================================================
# Scatter plots of true Cij against posterior median for all nine models
# (rows = prior, columns = sparsity), shown separately for the truly non-zero
# coefficients and for the true zeros.

# NOTE: `posterior_Cij` and `true_Cij` are rebuilt below, overwriting the
# versions created earlier in the script. The contents are the same files, so
# this repeat is redundant and the two blocks can be removed if the earlier
# versions are still in memory.

# Posterior median Cij matrices for the nine models
posterior_Cij <- list(
  uniform_0 = readRDS(
    "posterior_Cij_uniform_0pct_sparsity.rds"
  ),
  uniform_25 = readRDS(
    "posterior_Cij_uniform_25pct_sparsity.rds"
  ),
  uniform_50 = readRDS(
    "posterior_Cij_uniform_50pct_sparsity.rds"
  ),
  
  horseshoe_0 = readRDS(
    "posterior_Cij_horseshoe_0pct_sparsity.rds"
  ),
  horseshoe_25 = readRDS(
    "posterior_Cij_horseshoe_25pct_sparsity.rds"
  ),
  horseshoe_50 = readRDS(
    "posterior_Cij_horseshoe_50pct_sparsity.rds"
  ),
  
  lasso_0 = readRDS(
    "posterior_Cij_lasso_0pct_sparsity.rds"
  ),
  lasso_25 = readRDS(
    "posterior_Cij_lasso_25pct_sparsity.rds"
  ),
  lasso_50 = readRDS(
    "posterior_Cij_lasso_50pct_sparsity.rds"
  )
)

# True Cij matrices (repeat of the earlier definition)
true_Cij <- list(
  `0`  = readRDS("truth_Cij_0pct_sparsity.rds")[[1]],
  `25` = readRDS("truth_Cij_25pct_sparsity.rds")[[1]],
  `50` = readRDS("truth_Cij_50pct_sparsity.rds")[[1]]
)



priors <- c("uniform", "horseshoe", "lasso")
sparsities <- c(0, 25, 50)

# Build one long data frame with the true value and posterior median of every
# coefficient, for every prior x sparsity combination (nine models x 200
# coefficients)
Cij_recovery_all <- bind_rows(
  lapply(priors, function(p) {
    
    bind_rows(
      lapply(sparsities, function(s) {
        
        model_name <- paste0(p, "_", s)
        
        true_matrix <- true_Cij[[as.character(s)]]
        estimated_matrix <- posterior_Cij[[model_name]]
        
        # Safety check: both matrices must be 20 x 10
        stopifnot(
          identical(dim(true_matrix), c(20L, 10L)),
          identical(dim(estimated_matrix), c(20L, 10L))
        )
        
        data.frame(
          prior = p,
          sparsity = s,
          true_Cij = as.vector(true_matrix),
          posterior_median = as.vector(estimated_matrix)
        )
      })
    )
  })
) %>%
  mutate(
    # Label each coefficient by whether its TRUE value is exactly zero
    interaction_type = ifelse(
      true_Cij == 0,
      "True zero",
      "True non-zero"
    ),
    # Display labels and order for the facets
    prior = factor(
      prior,
      levels = priors,
      labels = c("Uniform", "Horseshoe", "LASSO")
    ),
    sparsity = factor(
      sparsity,
      levels = c(0, 25, 50),
      labels = c("0%", "25%", "50%")
    )
  )

# Check how many true zero/non-zero interactions
# are included in each panel
table(
  Cij_recovery_all$prior,
  Cij_recovery_all$sparsity,
  Cij_recovery_all$interaction_type
)


#-------------------------------------------------------------------------------
# True non-zero coefficients: true value vs posterior median
#-------------------------------------------------------------------------------
# Points on the black dashed 1:1 line are perfectly recovered. The red dashed
# line is a linear regression of estimate on truth in each panel: a slope below
# 1 indicates shrinkage of the estimates towards zero.

Cij_nonzero_plot <- ggplot(
  Cij_recovery_all %>%
    filter(interaction_type == "True non-zero"),
  aes(
    x = true_Cij,
    y = posterior_median,
    colour = prior
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  geom_point(
    size = 1.5,
    alpha = 0.7
  ) +
  # Fitted regression line for each panel
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "red",
    linetype = "dashed",
    linewidth = 1.1
  ) +
  
  scale_colour_manual(
    values = c(
      "Uniform" = "#0072B2",
      "Horseshoe" = "#E69F00",
      "LASSO" = "#009E73"
    )
  ) +
  facet_grid(prior ~ sparsity) +
  coord_equal(
    xlim = c(-1, 1),
    ylim = c(-1, 1)
  ) +
  scale_x_continuous(
    breaks = seq(-1, 1, 0.5)
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, 0.5)
  ) +
  labs(
    x = expression("True interaction coefficient (" * C[ij] * ")"),
    y = expression("Posterior median (" * C[ij] * ")")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    # Increase spacing between rows
    panel.spacing.y = grid::unit(1.2, "cm")
  )

print(Cij_nonzero_plot)





#-------------------------------------------------------------------------------
# True zero coefficients: posterior medians
#-------------------------------------------------------------------------------
# Every true value here is 0, so all points lie on x = 0 and the vertical
# spread shows how far the posterior medians are from zero. The dashed
# horizontal line at y = 0 marks perfect recovery. The 0% sparsity panels are
# empty because that scenario has no true zeros.

Cij_zero_plot <- ggplot(
  Cij_recovery_all %>%
    filter(interaction_type == "True zero"),
  aes(
    x = true_Cij,
    y = posterior_median,
    colour = prior
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  
  geom_point(
    shape = 17,
    size = 1.8,
    alpha = 0.7
  ) +
  
  scale_colour_manual(
    values = c(
      "Uniform" = "#0072B2",
      "Horseshoe" = "#E69F00",
      "LASSO" = "#009E73"
    )
  ) +
  
  facet_grid(prior ~ sparsity) +
  
  scale_x_continuous(
    limits = c(-0.1, 0.1),
    breaks = c(-0.1, 0, 0.1)
  ) +
  
  labs(
    x = expression("True interaction coefficient (" * C[ij] * ")"),
    y = expression("Posterior median (" * C[ij] * ")")
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(Cij_zero_plot)