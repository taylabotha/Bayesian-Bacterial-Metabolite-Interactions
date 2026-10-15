#===============================================================================
# Cij_recovery_analysis.R
#
# Purpose:
#   Compare how well three Bayesian priors (Uniform, Horseshoe, LASSO) recover
#   the true microbial interaction matrix (Cij) in simulated microbiota data.
#   Each prior is evaluated at three true sparsity levels (0%, 25%, 50%),
#   giving 9 fitted Stan models in total (10 species x 3 split, 1000 iterations
#   with 500 warm-up).
#
# Workflow:
#   1. Load the saved Stan fits and the true (simulated) Cij matrices
#   2. Extract posterior median Cij estimates from each fit
#   3. Compute accuracy of the estimates (Mean Absolute Error, MAE)
#   4. Classify interactions as present/absent and compute sensitivity,
#      specificity, accuracy and precision
#   5. Plot sensitivity and specificity by prior and sparsity
#      (line plot and bar plot versions)
#
# Inputs:  .rds files containing saved rstan fits and true Cij matrices
# Outputs: data frames `accuracy_results` and `recovery_results`, and the
#          sensitivity/specificity plots
#===============================================================================


#===============================================================================
# 0. SETUP
#===============================================================================

# Libraries needed
library(deSolve)    # ODE solvers (used when simulating the microbial dynamics)
library(rstan)      # Bayesian inference with Stan; used here to extract draws
library(tidyr)      # pivot_longer() for reshaping results
library(ggplot2)    # plotting
library(dplyr)      # data manipulation (bind_rows, select, mutate, %>%)
library(posterior)  # tools for working with posterior draws
library(stringr)    # string utilities

# Set the working directory to the folder containing the saved .rds files.
# NOTE: this path is specific to the author's machine; change it before running.
setwd("/Users/tayni/Downloads/thesis proj- microbiota")

# Install and load cowplot (only needed for combining plots).
# NOTE: install.packages() only needs to be run once; it can be commented out
# after the first run.
install.packages("cowplot")
library(cowplot)


#===============================================================================
# 1. LOAD SAVED FITS
#===============================================================================
# Read in all 9 saved Stan fit files (3 priors x 3 sparsity levels).
# File names encode the prior, the 10x3 design, the sparsity level and the
# sampler settings (1000 iterations - 500 warm-up).
# NOTE: some file names below are shortened with "......" as placeholders
# for the full file names.

# --- Uniform prior ---

# Uniform prior, 25% sparsity
uniform_25 <- readRDS("......25.rds")

# Original 10x3 Uniform 0% sparsity fit
uniform_0_correct <- readRDS(
  "C:/Users/tayni/Downloads/thesis proj- microbiota/uniform-0%spar-10-3split-1000-500/.......rds"
)

# Original 10x3 Uniform 50% sparsity fit
uniform_50_correct <- readRDS(
  "C:/Users/tayni/Downloads/thesis proj- microbiota/uniform-50%spar-10-3split-1000-500/......50.rds"
)

# --- Horseshoe prior (0%, 25%, 50% sparsity) ---
horseshoe_0  <- readRDS("horseshoe_fit_10x3_0%1000-500.rds")
horseshoe_25 <- readRDS("horseshoe_fit_10x3_25%1000-500.rds")
horseshoe_50 <- readRDS("horseshoe_fit_10x3_50%1000-500.rds")

# --- LASSO prior (0%, 25%, 50% sparsity) ---
lasso_0  <- readRDS("lasso_fit_10x3_0%1000-500.rds")
lasso_25 <- readRDS("lasso_fit_10x3_25%1000-500.rds")
lasso_50 <- readRDS("lasso_fit_10x3_50%1000-500.rds")

# Collect all 9 fits into a single named list.
# Each saved .rds object is a list whose first element ([[1]]) is the stanfit.
# Names follow the pattern <prior>_<sparsity>[_correct]; the regular
# expressions further down parse these names to recover prior and sparsity.
fits <- list(
  uniform_0_correct  = uniform_0_correct[[1]],
  uniform_25         = uniform_25[[1]],
  uniform_50_correct = uniform_50_correct[[1]],
  
  horseshoe_0        = horseshoe_0[[1]],
  horseshoe_25       = horseshoe_25[[1]],
  horseshoe_50       = horseshoe_50[[1]],
  
  lasso_0            = lasso_0[[1]],
  lasso_25           = lasso_25[[1]],
  lasso_50           = lasso_50[[1]]
)

# Load the true (simulated) interaction matrices used to generate the data,
# one for each sparsity level. These are the ground truth for evaluation.
truth_0  <- readRDS("truth_10x3_uni0%1000-500.rds")
truth_25 <- readRDS("truth_10x3_uni25%1000-500.rds")
truth_50 <- readRDS("truth_10x3_uni50%1000-500.rds")

# Store the true Cij matrices in a list indexed by sparsity level ("0", "25",
# "50") so they can be looked up by name for each model.
true_Cij <- list(
  "0"  = truth_0[[1]],
  "25" = truth_25[[1]],
  "50" = truth_50[[1]]
)


#===============================================================================
# 2. POSTERIOR MEDIAN Cij ESTIMATES
#===============================================================================

# Extract the posterior median Cij matrix directly from each fit.
# rstan::extract() returns a 3-D array of draws (draws x rows x columns);
# apply(..., c(2, 3), median) takes the median over draws for each matrix
# element, giving one point-estimate matrix per model.
posterior_Cij <- lapply(fits, function(fit) {
  
  posterior_samples <- rstan::extract(
    fit,
    pars = "Cij"
  )$Cij
  
  apply(
    posterior_samples,
    c(2, 3),
    median
  )
})


#===============================================================================
# 3. ACCURACY OF ESTIMATES (MAE)
#===============================================================================

# Calculate the Mean Absolute Error (MAE) between the true and estimated Cij
# for all nine models. Lower MAE = better recovery of interaction strengths.
accuracy_results <- lapply(names(fits), function(model_name) {
  
  # Parse the sparsity level (0, 25 or 50) from the model name
  sparsity <- as.numeric(
    sub(".*_(0|25|50)(_correct)?$", "\\1", model_name)
  )
  
  # Parse the prior name (uniform, horseshoe or lasso) from the model name
  prior <- sub(
    "_(0|25|50)(_correct)?$", "", model_name
  )
  
  # Flatten the true matrix for the matching sparsity level into a vector
  true_vals <- as.vector(
    true_Cij[[as.character(sparsity)]]
  )
  
  # Flatten the posterior median estimate matrix into a vector
  est_vals <- as.vector(
    posterior_Cij[[model_name]]
  )
  
  # One row of results per model
  data.frame(
    model = model_name,
    prior = prior,
    sparsity = sparsity,
    MAE = mean(abs(true_vals - est_vals))
  )
  
}) %>%
  bind_rows()   # stack the nine rows into one data frame

print(accuracy_results)


# NOTE: The two blocks above (posterior_Cij and accuracy_results) were repeated
# verbatim in the original script. The repeat produces identical output and is
# redundant; it is kept here, commented identically, so the original script is
# unchanged. It can safely be deleted.

# Extract posterior median Cij matrices directly from the fits (repeat)
posterior_Cij <- lapply(fits, function(fit) {
  
  posterior_samples <- rstan::extract(
    fit,
    pars = "Cij"
  )$Cij
  
  apply(
    posterior_samples,
    c(2, 3),
    median
  )
})

# Calculate MAE for all nine models (repeat)
accuracy_results <- lapply(names(fits), function(model_name) {
  
  sparsity <- as.numeric(
    sub(".*_(0|25|50)(_correct)?$", "\\1", model_name)
  )
  
  prior <- sub(
    "_(0|25|50)(_correct)?$", "", model_name
  )
  
  true_vals <- as.vector(
    true_Cij[[as.character(sparsity)]]
  )
  
  est_vals <- as.vector(
    posterior_Cij[[model_name]]
  )
  
  data.frame(
    model = model_name,
    prior = prior,
    sparsity = sparsity,
    MAE = mean(abs(true_vals - est_vals))
  )
  
}) %>%
  bind_rows()

print(accuracy_results)


#===============================================================================
# 4. INTERACTION RECOVERY (CLASSIFICATION METRICS)
#===============================================================================
# Treat recovery of the interaction network as a binary classification problem:
#   - "true positive"  = a real (non-zero) interaction that the model detects
#   - "true negative"  = a true zero that the model also estimates as ~zero

# Classification threshold: an estimated interaction is called "present" only
# if its absolute posterior median exceeds this value. Posterior medians are
# rarely exactly zero, so a threshold is needed to declare an interaction absent.
classification_threshold <- 0.1

# Calculate recovery metrics for all nine models
recovery_results <- lapply(names(fits), function(model_name) {
  
  # Extract sparsity from model name
  sparsity <- as.numeric(
    sub(".*_(0|25|50)(_correct)?$", "\\1", model_name)
  )
  
  # Extract prior from model name
  prior <- sub(
    "_(0|25|50)(_correct)?$", "", model_name
  )
  
  # True interaction coefficients
  true_vals <- as.vector(
    true_Cij[[as.character(sparsity)]]
  )
  
  # Posterior median estimates
  est_vals <- as.vector(
    posterior_Cij[[model_name]]
  )
  
  # Safety check: true and estimated matrices must be the same size
  stopifnot(length(true_vals) == length(est_vals))
  
  # Classify true and estimated interactions:
  # truth is "non-zero" if the simulated coefficient is not exactly zero;
  # the estimate is "non-zero" if it exceeds the classification threshold
  true_nonzero <- abs(true_vals) > 0
  
  estimated_nonzero <- abs(est_vals) > classification_threshold
  
  # Confusion matrix counts
  TP <- sum(true_nonzero & estimated_nonzero)     # real interaction, detected
  TN <- sum(!true_nonzero & !estimated_nonzero)   # true zero, estimated as zero
  FP <- sum(!true_nonzero & estimated_nonzero)    # true zero, wrongly detected
  FN <- sum(true_nonzero & !estimated_nonzero)    # real interaction, missed
  
  # Sensitivity (true positive rate): proportion of real interactions detected.
  # Returns NA if there are no true interactions.
  sensitivity <- if ((TP + FN) > 0) {
    TP / (TP + FN)
  } else {
    NA_real_
  }
  
  # Specificity (true negative rate): proportion of true zeros correctly
  # identified as zero. Returns NA if there are no true zeros.
  specificity <- if ((TN + FP) > 0) {
    TN / (TN + FP)
  } else {
    NA_real_
  }
  
  # Overall classification accuracy: proportion of all entries classified correctly
  accuracy <- (TP + TN) / length(true_vals)
  
  # Precision: proportion of detected interactions that are real.
  # Returns NA if nothing was detected.
  precision <- if ((TP + FP) > 0) {
    TP / (TP + FP)
  } else {
    NA_real_
  }
  
  # One row of results per model
  data.frame(
    model = model_name,
    prior = prior,
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

# Print results
print(recovery_results)


#===============================================================================
# 5. PLOTS: SENSITIVITY AND SPECIFICITY
#===============================================================================

#-------------------------------------------------------------------------------
# 5a. Line plot version
#-------------------------------------------------------------------------------

# Convert sensitivity and specificity into long format (one row per
# model-metric pair) so ggplot can facet by metric.
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
    # Standardise prior names and set a fixed display order / nice labels
    prior = factor(
      tolower(as.character(prior)),
      levels = c("uniform", "horseshoe", "lasso"),
      labels = c("Uniform", "Horseshoe", "LASSO")
    ),
    
    # Keep sparsity numeric (0, 25, 50) so it can be used on a continuous axis
    sparsity = as.numeric(
      gsub("%", "", as.character(sparsity))
    ),
    
    # Capitalise metric names for the panel titles
    metric = factor(
      metric,
      levels = c("sensitivity", "specificity"),
      labels = c("Sensitivity", "Specificity")
    )
  )

# Create two-panel line graph (Sensitivity | Specificity),
# with one coloured line per prior across the three sparsity levels
Sensitivity_Specificity_plot <- ggplot(
  sensitivity_specificity_df,
  aes(
    x = sparsity,
    y = value,
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
  
  # One panel per metric
  facet_wrap(
    ~ metric,
    ncol = 2
  ) +
  
  # Show the sparsity levels as percentages
  scale_x_continuous(
    breaks = c(0, 25, 50),
    labels = c("0%", "25%", "50%")
  ) +
  
  # Both metrics are proportions between 0 and 1
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  
  # Colour-blind-friendly palette (Okabe-Ito)
  scale_colour_manual(
    values = c(
      "Uniform" = "#0072B2",
      "Horseshoe" = "#D55E00",
      "LASSO" = "#009E73"
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


#-------------------------------------------------------------------------------
# 5b. Bar plot version
#-------------------------------------------------------------------------------
# NOTE: this block rebuilds the long-format data (sparsity is now a categorical
# factor instead of numeric) and OVERWRITES the objects
# `sensitivity_specificity_df` and `Sensitivity_Specificity_plot` from 5a.
# Save the line plot first (e.g. with ggsave()) if you want to keep both.

# Convert sensitivity and specificity into long format
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
    prior = factor(
      tolower(as.character(prior)),
      levels = c("uniform", "horseshoe", "lasso"),
      labels = c("Uniform", "Horseshoe", "LASSO")
    ),
    
    # Sparsity as an ordered categorical variable with "%" labels,
    # so it works as the x-axis for grouped bars
    sparsity = factor(
      paste0(gsub("%", "", as.character(sparsity)), "%"),
      levels = c("0%", "25%", "50%")
    ),
    
    metric = factor(
      metric,
      levels = c("sensitivity", "specificity"),
      labels = c("Sensitivity", "Specificity")
    )
  )

# Create two-panel bar graph: bars are grouped by sparsity level and
# coloured by prior, with one panel each for sensitivity and specificity
Sensitivity_Specificity_plot <- ggplot(
  sensitivity_specificity_df,
  aes(
    x = sparsity,
    y = value,
    fill = prior
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.8),   # place the three priors side by side
    width = 0.7,
    color = "black",                          # thin black outline on each bar
    linewidth = 0.3,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ metric,
    ncol = 2
  ) +
  
  # Same colour palette as the line plot, except Horseshoe is orange (#E69F00)
  scale_fill_manual(
    values = c(
      "Uniform" = "#0072B2",
      "Horseshoe" = "#E69F00",
      "LASSO" = "#009E73"
    )
  ) +
  
  # Axis runs slightly above 1 so bars at 1.0 are not clipped;
  # bars start flush with the x-axis (no expansion)
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = "Sparsity",
    y = "Proportion correctly classified",
    fill = "Prior"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      face = "bold",
      size = 11
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()    # remove vertical gridlines for cleaner bars
  )

print(Sensitivity_Specificity_plot)