#===============================================================================
# MCMC_diagnostics_and_posterior_plots.R
#
# Purpose:
#   Assess MCMC convergence and visualise the posterior distributions of the
#   interaction coefficients (Cij) for the nine fitted Stan models:
#   3 priors (Uniform, Horseshoe, LASSO) x 3 sparsity levels (0%, 25%, 50%),
#   all from the 10x3 simulated community (1000 iterations, 500 warm-up).
#
# Contents:
#   1. Load the saved Stan fits
#   2. Traceplots: 18 panels (3 priors x 3 sparsity levels x 2 coefficients)
#      -> Combined_18_Traceplots_3x6.png
#   3. Convergence diagnostics: Rhat and n_eff summaries, and an Rhat histogram
#   4. Posterior density plots with the true value marked: 18 panels
#      -> Posterior_distributions_10x3_3x6.png
#
# Inputs:  .rds files with saved rstan fits and the true Cij matrices
# Outputs: two high-resolution PNG figures, plus Rhat / n_eff summaries printed
#          to the console
#===============================================================================


#===============================================================================
# 0. SETUP
#===============================================================================

# Libraries needed
library(deSolve)    # ODE solvers (used when simulating the microbial dynamics)
library(rstan)      # Bayesian inference with Stan; used for extract/traceplot
library(tidyr)      # pivot_longer() for reshaping posterior draws
library(ggplot2)    # plotting
library(dplyr)      # data manipulation (bind_rows, filter, mutate, %>%)
library(posterior)  # tools for working with posterior draws
library(stringr)    # string utilities (str_detect)

# Set the working directory to the folder containing the saved .rds files.
# NOTE: this path is specific to the author's machine; change it before running.
setwd("/Users/tayni/Downloads/thesis proj- microbiota")

# Install and load cowplot (used to arrange multiple plots into one figure).
# NOTE: install.packages() only needs to be run once; it can be commented out
# after the first run.
install.packages("cowplot")
library(cowplot)


#===============================================================================
# 1. LOAD SAVED FITS
#===============================================================================

# Read in all 9 Stan fit files (3 priors x 3 sparsity levels).
# File names encode the prior, the 10x3 design, the sparsity level and the
# sampler settings (1000 iterations - 500 warm-up).
# NOTE: some file names are shortened with "......" as placeholders for the
# full names.

# --- Uniform prior (0%, 25%, 50% sparsity) ---
uniform_0  <- readRDS(".......rds")
uniform_25 <- readRDS("......25.rds")
uniform_50 <- readRDS("......50.rds")

# --- Horseshoe prior (0%, 25%, 50% sparsity) ---
horseshoe_0  <- readRDS("horseshoe_fit_10x3_0%1000-500.rds")
horseshoe_25 <- readRDS("horseshoe_fit_10x3_25%1000-500.rds")
horseshoe_50 <- readRDS("horseshoe_fit_10x3_50%1000-500.rds")

# --- LASSO prior (0%, 25%, 50% sparsity) ---
lasso_0  <- readRDS("lasso_fit_10x3_0%1000-500.rds")
lasso_25 <- readRDS("lasso_fit_10x3_25%1000-500.rds")
lasso_50 <- readRDS("lasso_fit_10x3_50%1000-500.rds")

# Collect all fits into one named list. Each saved .rds object is a list whose
# first element ([[1]]) is the stanfit object.
fits <- list(
  uniform_0      = uniform_0[[1]],
  uniform_25     = uniform_25[[1]],
  uniform_50     = uniform_50[[1]],
  
  horseshoe_0    = horseshoe_0[[1]],
  horseshoe_25   = horseshoe_25[[1]],
  horseshoe_50   = horseshoe_50[[1]],
  
  lasso_0        = lasso_0[[1]],
  lasso_25       = lasso_25[[1]],
  lasso_50       = lasso_50[[1]]
)

# Load the original 10x3 Uniform 0% fit
uniform_0_correct <- readRDS(
  "C:/Users/tayni/Downloads/thesis proj- microbiota/uniform-0%spar-10-3split-1000-500/.......rds"
)

# Load the original 10x3 Uniform 50% fit
uniform_50_correct <- readRDS(
  "C:/Users/tayni/Downloads/thesis proj- microbiota/uniform-50%spar-10-3split-1000-500/......50.rds"
)

# Replace the incorrect Uniform 0% and 50% models in `fits` with the original
# (correct) fits loaded above
fits[["uniform_0"]]  <- uniform_0_correct[[1]]
fits[["uniform_50"]] <- uniform_50_correct[[1]]


#===============================================================================
# 2. TRACEPLOTS
#===============================================================================
# Traceplots of the post-warm-up chains for two selected coefficients in each
# of the nine models (18 panels), arranged as a 3 x 6 grid:
#   rows    = prior (Uniform, Horseshoe, LASSO)
#   columns = sparsity level (0%, 25%, 50%) x coefficient (C[1,1], C[10,2])

# Libraries (already loaded above; repeated here so the section can be run alone)
library(rstan)
library(ggplot2)
library(cowplot)

#-------------------------------------------------------------------------------
# 2a. List all nine fitted models
#-------------------------------------------------------------------------------
# Named list with display labels. The order (Uniform, Horseshoe, LASSO; each at
# 0%, 25%, 50%) determines the order of panels in the final figure.

fits_trace <- list(
  "Uniform 0%"   = fits[["uniform_0"]],
  "Uniform 25%"  = fits[["uniform_25"]],
  "Uniform 50%"  = fits[["uniform_50"]],
  
  "Horseshoe 0%"  = fits[["horseshoe_0"]],
  "Horseshoe 25%" = fits[["horseshoe_25"]],
  "Horseshoe 50%" = fits[["horseshoe_50"]],
  
  "LASSO 0%"  = fits[["lasso_0"]],
  "LASSO 25%" = fits[["lasso_25"]],
  "LASSO 50%" = fits[["lasso_50"]]
)

#-------------------------------------------------------------------------------
# 2b. Select the same two coefficients for all models
#-------------------------------------------------------------------------------
# Cij[1,1] is plotted on a 0 to 1 scale and Cij[10,2] on a -1 to 0 scale
# (see make_traceplot() below, which assumes these signs).

selected_Cij <- c("Cij[1,1]", "Cij[10,2]")

#-------------------------------------------------------------------------------
# 2c. Function to create each traceplot
#-------------------------------------------------------------------------------
# Arguments:
#   fit        - a stanfit object
#   parameter  - the parameter name to plot (e.g. "Cij[1,1]")
#   model_name - display label for the model (used in the title)

make_traceplot <- function(fit, parameter, model_name) {
  
  # Set the y-axis range according to the coefficient's sign: Cij[1,1] is
  # positive in these simulations, Cij[10,2] is negative. A fixed range keeps
  # all panels for the same coefficient on a comparable scale.
  if (parameter == "Cij[1,1]") {
    y_limits <- c(0, 1)
  } else {
    y_limits <- c(-1, 0)
  }
  
  # rstan::traceplot() returns a ggplot object, so ggplot layers can be added.
  # inc_warmup = FALSE shows only the post-warm-up samples.
  traceplot(
    fit,
    pars = parameter,
    inc_warmup = FALSE
  ) +
    
    # Zoom the y-axis without dropping any data
    coord_cartesian(ylim = y_limits) +
    
    labs(
      title = paste(model_name, "-", parameter),
      x = "Iteration",
      y = parameter
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      axis.title = element_text(
        face = "bold",
        size = 9
      ),
      axis.text = element_text(size = 8),
      legend.position = "none"   # chain colour legend is omitted
    )
}

#-------------------------------------------------------------------------------
# 2d. Generate all 18 traceplots
#-------------------------------------------------------------------------------
# Loops over models (outer) and coefficients (inner), so the list order is:
# Uniform 0% [1,1], Uniform 0% [10,2], Uniform 25% [1,1], ... LASSO 50% [10,2].

traceplots <- list()

for (model_name in names(fits_trace)) {
  
  for (parameter in selected_Cij) {
    
    traceplots[[length(traceplots) + 1]] <- make_traceplot(
      fit = fits_trace[[model_name]],
      parameter = parameter,
      model_name = model_name
    )
    
  }
}

#-------------------------------------------------------------------------------
# 2e. Change subplot titles and remove individual y-axis labels
#-------------------------------------------------------------------------------
# The prior is labelled once per row in the final figure, so each panel title
# only needs the sparsity level and coefficient.

subplot_titles <- c(
  "0%: C[1,1]",
  "0%: C[10,2]",
  "25%: C[1,1]",
  "25%: C[10,2]",
  "50%: C[1,1]",
  "50%: C[10,2]"
)

for (i in seq_along(traceplots)) {
  
  traceplots[[i]] <- traceplots[[i]] +
    
    # ((i - 1) %% 6) + 1 cycles through the 6 titles, repeating for each prior
    labs(
      title = subplot_titles[((i - 1) %% 6) + 1],
      x = "Iteration",
      y = NULL
    ) +
    
    # Tick marks every 200 iterations (post-warm-up chain length is 500)
    scale_x_continuous(
      breaks = seq(0, 1000, by = 200)
    ) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 9
      ),
      
      # Remove vertical Cij label
      axis.title.y = element_blank(),
      
      # Make axis text smaller to prevent overlap
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      
      axis.title.x = element_text(
        size = 8,
        face = "bold"
      ),
      
      # Reduce unnecessary space around each plot
      plot.margin = margin(3, 1, 3, 1)
    )
}

#-------------------------------------------------------------------------------
# 2f. Combine traceplots: 3 rows x 6 columns
#-------------------------------------------------------------------------------

library(cowplot)

# Each row holds the six panels for one prior
row1 <- plot_grid(
  plotlist = traceplots[1:6],    # Uniform
  ncol = 6
)

row2 <- plot_grid(
  plotlist = traceplots[7:12],   # Horseshoe
  ncol = 6
)

row3 <- plot_grid(
  plotlist = traceplots[13:18],  # LASSO
  ncol = 6
)

# Helper: add a rotated prior label on the left-hand side of a row of plots.
# Draws the label in a narrow strip at the left and the plot in the remaining
# width.
add_prior_label <- function(plot, prior_name) {
  
  ggdraw() +
    
    draw_label(
      prior_name,
      x = 0.015,
      y = 0.5,
      angle = 90,
      fontface = "bold",
      size = 12
    ) +
    
    draw_plot(
      plot,
      x = 0.035,
      y = 0,
      width = 0.965,
      height = 1
    )
}

row1 <- add_prior_label(row1, "Uniform")
row2 <- add_prior_label(row2, "Horseshoe")
row3 <- add_prior_label(row3, "LASSO")

# Stack the three labelled rows into the final figure
combined_traceplots <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1
)

print(combined_traceplots)

#-------------------------------------------------------------------------------
# 2g. Save high-resolution figure
#-------------------------------------------------------------------------------

ggsave(
  filename = "Combined_18_Traceplots_3x6.png",
  plot = combined_traceplots,
  width = 24,
  height = 10,
  dpi = 300,
  limitsize = FALSE
)


#===============================================================================
# 3. CONVERGENCE DIAGNOSTICS (Rhat AND n_eff)
#===============================================================================

# NOTE: This block reloads the nine fits and rebuilds `fits`. It repeats
# Section 1, EXCEPT that the Uniform 0% and 50% fits are NOT replaced with
# `uniform_0_correct` / `uniform_50_correct` afterwards. If the first-pass
# Uniform 0% / 50% files are the "incorrect" models mentioned in Section 1,
# the diagnostics and posterior plots below will be based on them.
# Consider deleting this block (the `fits` list from Section 1 is still in
# memory) or re-applying the replacement here.

# Read in all 9 Stan fit files
uniform_0  <- readRDS(".......rds")
uniform_25 <- readRDS("......25.rds")
uniform_50 <- readRDS("......50.rds")

horseshoe_0  <- readRDS("horseshoe_fit_10x3_0%1000-500.rds")
horseshoe_25 <- readRDS("horseshoe_fit_10x3_25%1000-500.rds")
horseshoe_50 <- readRDS("horseshoe_fit_10x3_50%1000-500.rds")

lasso_0  <- readRDS("lasso_fit_10x3_0%1000-500.rds")
lasso_25 <- readRDS("lasso_fit_10x3_25%1000-500.rds")
lasso_50 <- readRDS("lasso_fit_10x3_50%1000-500.rds")

fits <- list(
  uniform_0      = uniform_0[[1]],
  uniform_25     = uniform_25[[1]],
  uniform_50     = uniform_50[[1]],
  
  horseshoe_0    = horseshoe_0[[1]],
  horseshoe_25   = horseshoe_25[[1]],
  horseshoe_50   = horseshoe_50[[1]],
  
  lasso_0        = lasso_0[[1]],
  lasso_25       = lasso_25[[1]],
  lasso_50       = lasso_50[[1]]
)

# Update the Uniform 0% and 50% entries in the traceplot list from `fits`
fits_trace[["Uniform 0%"]]  <- fits[["uniform_0"]]
fits_trace[["Uniform 50%"]] <- fits[["uniform_50"]]

# Extract Rhat and n_eff for every parameter in every model.
# rstan::summary(f)$summary is a matrix with one row per parameter; the
# "Rhat" and "n_eff" columns hold the convergence diagnostics.
# bind_rows(.id = "model") stacks the nine data frames and records the model
# name (the name of each element of `fits`) in a new "model" column.
diagnostics_df <- lapply(fits, function(f) {
  
  summ <- rstan::summary(f)$summary
  
  data.frame(
    parameter = rownames(summ),
    Rhat = summ[, "Rhat"],
    neff = summ[, "n_eff"]
  )
  
}) %>%
  bind_rows(.id = "model")

# Keep only the interaction coefficients (Cij[i,j]) and the sigma_ar parameter
selected_diagnostics <- diagnostics_df %>%
  filter(
    str_detect(parameter, "Cij\\[") |
      parameter == "sigma_ar"
  )

#-------------------------------------------------------------------------------
# 3a. Print overall minimum, maximum and range
#-------------------------------------------------------------------------------
# Summaries pool the selected parameters across all nine models.
# Rhat close to 1 (below 1.01 here) indicates convergence; a larger n_eff
# indicates more independent information in the chains.

cat("\n==============================")
cat("\nRHAT RESULTS - ALL 9 MODELS")
cat("\n==============================\n")

cat("Minimum Rhat:",
    min(selected_diagnostics$Rhat, na.rm = TRUE), "\n")

cat("Maximum Rhat:",
    max(selected_diagnostics$Rhat, na.rm = TRUE), "\n")

cat("Range of Rhat:",
    diff(range(selected_diagnostics$Rhat, na.rm = TRUE)), "\n")


cat("\n==============================")
cat("\nN_EFF RESULTS - ALL 9 MODELS")
cat("\n==============================\n")

cat("Minimum n_eff:",
    min(selected_diagnostics$neff, na.rm = TRUE), "\n")

cat("Maximum n_eff:",
    max(selected_diagnostics$neff, na.rm = TRUE), "\n")

cat("Range of n_eff:",
    diff(range(selected_diagnostics$neff, na.rm = TRUE)), "\n")

#-------------------------------------------------------------------------------
# 3b. Rhat histogram (all models, coloured by prior)
#-------------------------------------------------------------------------------
# NOTE: The original section header said "9 PANEL FIGURE", but the code below
# produces a single histogram with the three priors overlaid.

library(ggplot2)
library(dplyr)

# Extract the prior from the model names. Works with either "uniform_0" or
# "Uniform - 0%" style names (matching is case-insensitive).
selected_diagnostics <- selected_diagnostics %>%
  mutate(
    prior = case_when(
      grepl("uniform", as.character(model), ignore.case = TRUE) ~ "Uniform",
      grepl("horseshoe", as.character(model), ignore.case = TRUE) ~ "Horseshoe",
      grepl("lasso", as.character(model), ignore.case = TRUE) ~ "LASSO"
    ),
    # Fix the display order of the priors
    prior = factor(prior, levels = c("Uniform", "Horseshoe", "LASSO"))
  )

# Combined Rhat histogram: one overlaid, semi-transparent histogram per prior
Rhat_plot <- ggplot(
  selected_diagnostics %>% filter(is.finite(Rhat)),   # drop NA / Inf values
  aes(x = Rhat, fill = prior)
) +
  
  geom_histogram(
    binwidth = 0.0005,
    position = "identity",   # overlay rather than stack
    alpha = 0.5,
    color = "black"
  ) +
  
  # Rhat convergence threshold (values below 1.01 indicate convergence)
  geom_vline(
    xintercept = 1.01,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Colour-blind-friendly palette
  scale_fill_manual(
    values = c(
      "Uniform" = "#0072B2",
      "Horseshoe" = "#E69F00",
      "LASSO" = "#009E73"
    )
  ) +
  
  labs(
    x = "Rhat",
    y = "Number of Parameters",
    fill = "Prior"
  ) +
  
  theme_minimal() +
  
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

print(Rhat_plot)

#-------------------------------------------------------------------------------
# 3c. n_eff figure
#-------------------------------------------------------------------------------
# NOTE: The original script has an "N_EFF - 9 PANEL FIGURE" header here but no
# code beneath it. The n_eff values are computed in `selected_diagnostics`
# ($neff) and summarised above; a plot can be added here if needed.


#===============================================================================
# 4. POSTERIOR DISTRIBUTIONS - SMALL COMMUNITY (10 x 3)
#===============================================================================
# Posterior density plots for the same two coefficients, with the true
# (simulated) value marked by a red dashed line:
# 3 priors x 3 sparsity levels x 2 coefficients = 18 panels.

#-------------------------------------------------------------------------------
# 4a. File names of all nine fits
#-------------------------------------------------------------------------------
# NOTE: `fit_files` is defined here but never used in this script (the fits are
# loaded individually above). Kept as a reference for the file naming scheme.

fit_files <- c(
  "uniform_0"    = "uniform_fit_10x3_0%1000-500.rds",
  "horseshoe_0"  = "horseshoe_fit_10x3_0%1000-500.rds",
  "lasso_0"      = "lasso_fit_10x3_0%1000-500.rds",
  
  "uniform_25"   = "uniform_fit_10x3_25%1000-500.rds",
  "horseshoe_25" = "horseshoe_fit_10x3_25%1000-500.rds",
  "lasso_25"     = "lasso_fit_10x3_25%1000-500.rds",
  
  "uniform_50"   = "uniform_fit_10x3_50%1000-500.rds",
  "horseshoe_50" = "horseshoe_fit_10x3_50%1000-500.rds",
  "lasso_50"     = "lasso_fit_10x3_50%1000-500.rds"
)

#-------------------------------------------------------------------------------
# 4b. Load the true Cij matrices
#-------------------------------------------------------------------------------
# One true interaction matrix per sparsity level, used to mark the true value
# on each posterior plot.

truth_0  <- readRDS("truth_10x3_uni0%1000-500.rds")
truth_25 <- readRDS("truth_10x3_uni25%1000-500.rds")
truth_50 <- readRDS("truth_10x3_uni50%1000-500.rds")

# Check what each file contains
str(truth_0)
str(truth_25)
str(truth_50)

# Store the true matrices in a list indexed by sparsity level ("0", "25", "50")
true_Cij <- list(
  "0"  = truth_0[[1]],
  "25" = truth_25[[1]],
  "50" = truth_50[[1]]
)

#-------------------------------------------------------------------------------
# 4c. Extract posterior draws (first version)
#-------------------------------------------------------------------------------
# NOTE: `get_posterior()` is defined twice in this script. This first version is
# superseded by the second definition in Section 4d, which overwrites it, so
# only the second version is used for the plots. This block is kept for
# reference and as a quick check of the number of draws per panel.
#
# For a model name like "horseshoe_25", this function:
#   - parses the prior ("horseshoe") and sparsity ("25") from the name
#   - extracts all posterior draws of the Cij matrix
#   - returns a long data frame of draws for Cij[1,1] and Cij[10,2], with the
#     true value, prior label and sparsity level attached

get_posterior <- function(model_name) {
  
  fit <- fits[[model_name]]
  
  # Split the model name into prior and sparsity (e.g. "lasso_50")
  model_parts <- strsplit(model_name, "_")[[1]]
  
  prior_name <- model_parts[1]
  sparsity <- model_parts[2]
  
  # Convert prior name to the display label before creating the data frame
  prior_label <- switch(
    tolower(prior_name),
    "uniform" = "Uniform",
    "horseshoe" = "Horseshoe",
    "lasso" = "LASSO"
  )
  
  if (is.null(prior_label)) {
    stop(paste("Unknown prior:", prior_name))
  }
  
  # True Cij matrix for this sparsity level
  true_matrix <- true_Cij[[sparsity]]
  
  # Posterior draws of the whole Cij matrix: a (draws x rows x columns) array
  Cij_draws <- rstan::extract(
    fit,
    pars = "Cij",
    permuted = TRUE
  )$Cij
  
  n_draws <- dim(Cij_draws)[1]
  
  # Long data frame: all draws for Cij[1,1], followed by all draws for Cij[10,2]
  posterior_df <- data.frame(
    parameter = rep(
      c("Cij[1,1]", "Cij[10,2]"),
      each = n_draws
    ),
    
    value = c(
      Cij_draws[, 1, 1],
      Cij_draws[, 10, 2]
    ),
    
    true_value = rep(
      c(true_matrix[1, 1], true_matrix[10, 2]),
      each = n_draws
    ),
    
    prior = prior_label,
    sparsity = as.integer(sparsity)
  )
  
  return(posterior_df)
}

# Apply to all nine models and stack the results
posterior_selected <- dplyr::bind_rows(
  lapply(names(fits), get_posterior)
)

# Sanity check: number of draws for each prior / sparsity / coefficient
posterior_selected %>%
  dplyr::count(prior, sparsity, parameter)

#-------------------------------------------------------------------------------
# 4d. Extract posterior draws (second version, the one used for the plots)
#-------------------------------------------------------------------------------

# Coefficients to plot
selected_parameters <- c("Cij[1,1]", "Cij[10,2]")

# This version extracts the selected parameters directly, reshapes them to long
# format and looks up each coefficient's true value from its row/column index
# parsed from the parameter name (e.g. "Cij[10,2]" -> row 10, column 2).
get_posterior <- function(model_name) {
  
  fit <- fits[[model_name]]
  
  # Split the model name into prior and sparsity (e.g. "lasso_50")
  model_parts <- strsplit(model_name, "_")[[1]]
  
  prior <- model_parts[1]
  sparsity <- model_parts[2]
  
  # True Cij matrix for this sparsity level
  true_matrix <- true_Cij[[sparsity]]
  
  # Posterior draws of the selected coefficients, one column per parameter
  draws <- as.data.frame(
    rstan::extract(
      fit,
      pars = selected_parameters,
      permuted = TRUE
    )
  )
  
  posterior_df <- draws %>%
    # One row per draw per parameter
    tidyr::pivot_longer(
      cols = all_of(selected_parameters),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      # Row index i from "Cij[i,j]"
      row_idx = as.integer(
        sub("Cij\\[([0-9]+),([0-9]+)\\]", "\\1", parameter)
      ),
      
      # Column index j from "Cij[i,j]"
      col_idx = as.integer(
        sub("Cij\\[([0-9]+),([0-9]+)\\]", "\\2", parameter)
      ),
      
      # Look up the true value for each coefficient in the true matrix
      true_value = mapply(
        function(i, j) true_matrix[i, j],
        row_idx,
        col_idx
      ),
      
      # Convert the prior name to its display label
      prior = case_when(
        prior == "uniform" ~ "Uniform",
        prior == "horseshoe" ~ "Horseshoe",
        prior == "lasso" ~ "LASSO"
      ),
      
      # Sparsity level as a number (0, 25 or 50)
      sparsity = as.integer(sparsity)
    )
  
  posterior_df
}

# Apply to all nine models and stack the results
posterior_selected <- bind_rows(
  lapply(names(fits), get_posterior)
)

#-------------------------------------------------------------------------------
# 4e. Create the 18 individual posterior plots
#-------------------------------------------------------------------------------
# Loop order (prior -> sparsity -> coefficient) gives six panels per prior,
# matching the 3 x 6 layout below: 0% C[1,1], 0% C[10,2], 25% C[1,1], ...

posterior_plots <- list()

prior_order <- c("Uniform", "Horseshoe", "LASSO")
sparsity_order <- c(0, 25, 50)

for (prior_name in prior_order) {
  
  for (s in sparsity_order) {
    
    for (parameter_name in selected_parameters) {
      
      # Posterior draws for this prior / sparsity / coefficient
      plot_data <- posterior_selected %>%
        filter(
          prior == prior_name,
          sparsity == s,
          parameter == parameter_name
        )
      
      # True value of this coefficient (constant within the subset)
      true_val <- unique(plot_data$true_value)
      
      # Shorter label for the panel title ("Cij[1,1]" -> "C[1,1]")
      short_parameter <- gsub("Cij", "C", parameter_name)
      
      p <- ggplot(plot_data, aes(x = value)) +
        
        # Posterior density
        geom_density(
          fill = "skyblue",
          alpha = 0.5,
          color = "black",
          linewidth = 0.5
        ) +
        
        # True value marked with a red dashed line
        geom_vline(
          xintercept = true_val,
          color = "red",
          linetype = "dashed",
          linewidth = 0.8
        ) +
        
        # Axis titles are removed from individual panels; a shared x-axis title
        # is added to the combined figure
        labs(
          title = paste0(s, "%: ", short_parameter),
          x = NULL,
          y = NULL
        ) +
        
        theme_minimal() +
        
        theme(
          plot.title = element_text(
            hjust = 0.5,
            face = "bold",
            size = 9
          ),
          
          axis.title.x = element_text(
            face = "bold",
            size = 8
          ),
          
          axis.title.y = element_blank(),
          
          axis.text = element_text(size = 7),
          
          plot.margin = margin(3, 1, 3, 1)
        )
      
      posterior_plots[[length(posterior_plots) + 1]] <- p
    }
  }
}

#-------------------------------------------------------------------------------
# 4f. Combine posterior plots: 3 rows x 6 columns
#-------------------------------------------------------------------------------

library(cowplot)

# Each row holds the six panels for one prior
row1 <- plot_grid(
  plotlist = posterior_plots[1:6],    # Uniform
  ncol = 6
)

row2 <- plot_grid(
  plotlist = posterior_plots[7:12],   # Horseshoe
  ncol = 6
)

row3 <- plot_grid(
  plotlist = posterior_plots[13:18],  # LASSO
  ncol = 6
)

# Add prior labels vertically on the left
# (same helper as in the traceplot section, redefined here so this section can
# be run on its own)
add_prior_label <- function(plot, prior_name) {
  
  ggdraw() +
    
    draw_label(
      prior_name,
      x = 0.015,
      y = 0.5,
      angle = 90,
      fontface = "bold",
      size = 12
    ) +
    
    draw_plot(
      plot,
      x = 0.035,
      y = 0,
      width = 0.965,
      height = 1
    )
}

row1 <- add_prior_label(row1, "Uniform")
row2 <- add_prior_label(row2, "Horseshoe")
row3 <- add_prior_label(row3, "LASSO")

# Combine all three rows
combined_rows <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1
)

# Add one shared x-axis title at the bottom: the combined plot occupies the
# upper 94% of the canvas, leaving a strip at the bottom for the label
Posterior_plot <- ggdraw() +
  
  draw_plot(
    combined_rows,
    x = 0,
    y = 0.06,
    width = 1,
    height = 0.94
  ) +
  
  draw_label(
    "Interaction coefficient",
    x = 0.5,
    y = 0.025,
    fontface = "bold",
    size = 12
  )

print(Posterior_plot)

#-------------------------------------------------------------------------------
# 4g. Save high-resolution figure
#-------------------------------------------------------------------------------

ggsave(
  filename = "Posterior_distributions_10x3_3x6.png",
  plot = Posterior_plot,
  width = 24,
  height = 10,
  dpi = 300,
  bg = "white"    # white background (the canvas would otherwise be transparent)
)
