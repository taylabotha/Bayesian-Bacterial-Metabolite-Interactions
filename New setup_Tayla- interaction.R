#===============================================================================
# Cij SBC pipeline — compare uniform / horseshoe / lasso priors
# across 0%, 25%, 50% sparsity, with each prior fit to the SAME simulated
# data within a given sparsity level (only the prior differs, not the data).
#===============================================================================

library(deSolve)
library(rstan)
library(tidyr)
library(ggplot2)
library(dplyr)
library(posterior)
library(stringr)

options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1)

#--------------------------------------------------------------------------------------
# 1. ODE system used to generate the "true" simulated data (deSolve) — unchanged logic
#--------------------------------------------------------------------------------------
bacteria_metabolite_ode <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N_met <- length(S_t0)
    N_bac <- length(y0_bacteria)
    B_total <- sum(state[(N_met + 1):(N_met + N_bac)])
    
    dydt <- numeric(N_met + N_bac)
    
    for (i in 1:N_met) {
      dydt[i] <- 0
      for (j in 1:N_bac) {
        dydt[i] <- dydt[i] - Cij[i, j] * state[i] * state[N_met + j]
      }
    }
    
    for (j in 1:N_bac) {
      growth <- -d_bac
      for (i in 1:N_met) {
        growth <- growth + C0 * Cij[i, j] * state[i]
      }
      dydt[N_met + j] <- state[N_met + j] * growth * (k - B_total)
    }
    
    list(dydt)
  })
}

#--------------------------------------------------------------------------------------
# 2. Constants shared across every sparsity level and every prior
#--------------------------------------------------------------------------------------
N_met <- 20
N_bac <- 10
Reps  <- 6
Obs   <- 4
k     <- 14
delta <- 1e-6
C0    <- 0.1
d_bac <- 0.1
n_sims <- 1
t0 <- 0
tmax <- 4
ts_sim <- seq(1e-6, tmax, length.out = Obs)

add_noise <- FALSE
Sigma_met   <- runif(N_met, min = 0.1, max = 0.5)
Sigma_bac   <- runif(N_bac, min = 0.01, max = 0.05)
Sigma_noise <- c(Sigma_met, Sigma_bac)

consumed_metabolites <- c(1:14)
produced_metabolites <- c(15:20)

# The three sparsity conditions
sparsity_levels <- c(0, 0.25, 0.50)

# The three priors to compare (must match names(Stan_models) below)
prior_names <- c("uniform", "horseshoe", "lasso")

set.seed(42)

#--------------------------------------------------------------------------------------
# 3. Helper: simulate ONE Cij matrix (at a given sparsity) + ODE-generated data.
#    Called once per sparsity level, then reused for all three priors at that
#    sparsity level, so the comparison isolates the effect of the prior.
#--------------------------------------------------------------------------------------
simulate_condition <- function(proportion_sparsity) {
  
  Cij_list             <- vector("list", n_sims)
  S_t0_list_all         <- vector("list", n_sims)
  y0_bacteria_list_all  <- vector("list", n_sims)
  y_obs_all <- array(0, dim = c(n_sims, Reps, Obs, N_met + N_bac))
  
  for (sim in 1:n_sims) {
    
    Cij_consumed <- matrix(runif(length(consumed_metabolites) * N_bac, 0, 1),
                           nrow = length(consumed_metabolites), ncol = N_bac)
    Cij_produced <- matrix(runif(length(produced_metabolites) * N_bac, -1, 0),
                           nrow = length(produced_metabolites), ncol = N_bac)
    Cij <- matrix(0, N_met, N_bac)
    Cij[consumed_metabolites, ] <- Cij_consumed
    Cij[produced_metabolites, ] <- Cij_produced
    
    n_zero <- round(proportion_sparsity * length(Cij))
    if (n_zero > 0) {
      zero_locations <- sample(length(Cij), n_zero)
      Cij[zero_locations] <- 0
    }
    
    Cij_list[[sim]] <- Cij
    
    S_t0_list <- replicate(Reps, {
      c(runif(length(consumed_metabolites), min = 0.5, max = 1),
        runif(length(produced_metabolites), min = 1e-3, max = 1e-2))
    }, simplify = FALSE)
    
    y0_bacteria_list <- replicate(Reps,
                                  runif(N_bac, min = 0.05, max = 0.2),
                                  simplify = FALSE)
    
    S_t0_list_all[[sim]]        <- S_t0_list
    y0_bacteria_list_all[[sim]] <- y0_bacteria_list
    
    for (r in 1:Reps) {
      state_init <- c(S_t0_list[[r]], y0_bacteria_list[[r]])
      
      parameters <- list(
        Cij = Cij,
        k = k,
        S_t0 = S_t0_list[[r]],
        y0_bacteria = y0_bacteria_list[[r]],
        C0 = C0,
        consumed_metabolites = consumed_metabolites,
        produced_metabolites = produced_metabolites,
        d_bac = d_bac
      )
      
      out <- ode(y = state_init, times = ts_sim, func = bacteria_metabolite_ode,
                 parms = parameters, method = "ode45")
      
      for (i in 1:(N_met + N_bac)) {
        true_vals <- out[, i + 1]
        
        if (!add_noise) {
          y_obs_all[sim, r, , i] <- true_vals + delta
        } else {
          noise <- rnorm(1, mean = 0, sd = Sigma_noise[i])
          y_obs_all[sim, r, , i] <- exp(log(true_vals + delta) + noise)
        }
      }
    }
  }
  
  list(Cij_list = Cij_list,
       S_t0_list_all = S_t0_list_all,
       y0_bacteria_list_all = y0_bacteria_list_all,
       y_obs_all = y_obs_all)
}

#--------------------------------------------------------------------------------------
# 4. Shared functions{} block (vectorized ODE — see earlier speed-up discussion:
#    replaces the 200-term nested loops with matrix algebra, and each model uses
#    ode_rk45_tol with loosened tolerances + a single reused ODE solve rather than
#    solving twice per iteration). All three models below use this identical block.
#--------------------------------------------------------------------------------------
functions_block <- "
functions {
  vector Bac_Met(real t, vector y, matrix Cij, real k, vector S_t0,
                 int[] consumed_metabolites, real C0) {
    int N_met = rows(Cij);
    int N_bac = cols(Cij);
    vector[N_met] y_met = y[1:N_met];
    vector[N_bac] y_bac = y[(N_met + 1):(N_met + N_bac)];
    real B_total = sum(y_bac);

    vector[N_met] dydt_met = -y_met .* (Cij * y_bac);

    vector[N_bac] growth = C0 * (to_row_vector(y_met[consumed_metabolites])
                                  * Cij[consumed_metabolites, ])' - 0.1;
    vector[N_bac] dydt_bac = y_bac .* growth * (k - B_total);

    return append_row(dydt_met, dydt_bac);
  }
}
"

data_block <- "
data {
  int<lower=1> R;
  int<lower=1> T;
  int<lower=1> N_met;
  int<lower=1> N_bac;
  array[T] real ts;
  array[R, T] vector[N_met + N_bac] y_obs;
  vector[N_met] S_t0[R];
  vector[N_bac] y0_bacteria[R];
  real t0;
  real<lower=0> k;
  int<lower=0> N_produced;
  int<lower=0> N_consumed;
  int produced_metabolites[N_produced];
  int consumed_metabolites[N_consumed];
  int<lower=0, upper=1> use_ar1;
  real<lower=0> delta;
  real<lower=0, upper=1> C0;
}
"

# The mu-solving block shared by all three models' transformed parameters{}
mu_block <- "
  array[R, T] vector[N_met + N_bac] mu;
  for (r in 1:R) {
    vector[N_met + N_bac] y0;
    y0[1:N_met] = S_t0[r];
    y0[(N_met + 1):(N_met + N_bac)] = y0_bacteria[r];
    mu[r] = ode_rk45_tol(Bac_Met, y0, t0, ts, 1e-4, 1e-4, 2000,
                          Cij, k, S_t0[r], consumed_metabolites, C0);
  }
"

# The likelihood loop shared by all three models' model{} block
likelihood_block <- "
  for (r in 1:R) {
    for (l in 1:(N_met + N_bac)) {
      vector[T] eps;
      for (t in 1:T) {
        eps[t] = log(y_obs[r, t, l] + delta) - log(mu[r, t, l] + delta);
      }
      eps ~ normal(0, sigma_ar);
    }
  }
"

# generated quantities shared by all three models (mu reused, not re-solved)
gq_block <- "
generated quantities {
  array[R, T] vector[N_met + N_bac] y_sim = mu;
  array[R, T, N_met + N_bac] real log_lik;

  for (r in 1:R) {
    for (t in 1:T) {
      for (l in 1:(N_met + N_bac)) {
        real eps = log(y_obs[r, t, l] + delta) - log(y_sim[r, t, l] + delta);
        log_lik[r, t, l] = normal_lpdf(eps | 0, sigma_ar);
      }
    }
  }
}
"

#--------------------------------------------------------------------------------------
# 5a. Uniform prior model
#--------------------------------------------------------------------------------------
Stan_code_uniform <- paste0(functions_block, data_block, "
parameters {
  matrix<lower=0, upper=1>[N_consumed, N_bac] Cij_pos;
  matrix<lower=-1, upper=0>[N_produced, N_bac] Cij_neg;
  real<lower=0> sigma_ar;
}

transformed parameters {
  matrix[N_met, N_bac] Cij;
  for (i in 1:N_consumed) {
    for (j in 1:N_bac) {
      Cij[i, j] = Cij_pos[i, j];
    }
  }
  for (i in 1:N_produced) {
    int idx = N_consumed + i;
    for (j in 1:N_bac) {
      Cij[idx, j] = Cij_neg[i, j];
    }
  }
", mu_block, "
}

model {
  sigma_ar ~ normal(0, 0.2);
", likelihood_block, "
}
", gq_block)

#--------------------------------------------------------------------------------------
# 5b. Horseshoe prior model
#--------------------------------------------------------------------------------------
Stan_code_horseshoe <- paste0(functions_block, data_block, "
parameters {
  real<lower=0> tau;
  matrix<lower=0>[N_consumed, N_bac] lambda_pos;
  matrix<lower=0>[N_produced, N_bac] lambda_neg;
  matrix<lower=0, upper=1>[N_consumed, N_bac] Cij_pos;
  matrix<lower=0, upper=1>[N_produced, N_bac] Cij_neg_mag;
  real<lower=0> sigma_ar;
}

transformed parameters {
  matrix[N_met, N_bac] Cij;
  for (i in 1:N_consumed) {
    for (j in 1:N_bac) {
      Cij[i, j] = Cij_pos[i, j];
    }
  }
  for (i in 1:N_produced) {
    int idx = N_consumed + i;
    for (j in 1:N_bac) {
      Cij[idx, j] = -Cij_neg_mag[i, j];
    }
  }
", mu_block, "
}

model {
  tau ~ cauchy(0, 1);
  to_vector(lambda_pos) ~ cauchy(0, 1);
  to_vector(lambda_neg) ~ cauchy(0, 1);

  for (i in 1:N_consumed) {
    for (j in 1:N_bac) {
      Cij_pos[i, j] ~ normal(0, tau * lambda_pos[i, j]) T[0, 1];
    }
  }
  for (i in 1:N_produced) {
    for (j in 1:N_bac) {
      Cij_neg_mag[i, j] ~ normal(0, tau * lambda_neg[i, j]) T[0, 1];
    }
  }

  sigma_ar ~ normal(0, 0.2);
", likelihood_block, "
}
", gq_block)

#--------------------------------------------------------------------------------------
# 5c. Lasso prior model
#--------------------------------------------------------------------------------------
Stan_code_lasso <- paste0(functions_block, data_block, "
parameters {
  matrix<lower=0, upper=1>[N_consumed, N_bac] Cij_pos;
  matrix<lower=0, upper=1>[N_produced, N_bac] Cij_neg_mag;
  real<lower=0> b_lasso;
  real<lower=0> sigma_ar;
}

transformed parameters {
  matrix[N_met, N_bac] Cij;
  for (i in 1:N_consumed) {
    for (j in 1:N_bac) {
      Cij[i, j] = Cij_pos[i, j];
    }
  }
  for (i in 1:N_produced) {
    int idx = N_consumed + i;
    for (j in 1:N_bac) {
      Cij[idx, j] = -Cij_neg_mag[i, j];
    }
  }
", mu_block, "
}

model {
  b_lasso ~ exponential(5);

  for (i in 1:N_consumed) {
    for (j in 1:N_bac) {
      Cij_pos[i, j] ~ double_exponential(0, b_lasso) T[0, 1];
    }
  }
  for (i in 1:N_produced) {
    for (j in 1:N_bac) {
      Cij_neg_mag[i, j] ~ double_exponential(0, b_lasso) T[0, 1];
    }
  }

  sigma_ar ~ normal(0, 0.2);
", likelihood_block, "
}
", gq_block)

#--------------------------------------------------------------------------------------
# 6. Compile each model ONCE (3 compiles total, not 9)
#--------------------------------------------------------------------------------------
Stan_models <- list(
  uniform   = stan_model(model_code = Stan_code_uniform),
  horseshoe = stan_model(model_code = Stan_code_horseshoe),
  lasso     = stan_model(model_code = Stan_code_lasso)
)

#--------------------------------------------------------------------------------------
# 7. Main loop: for each sparsity level, simulate data ONCE, then fit all three
#    priors to that SAME data. 3 sparsity levels x 3 priors = 9 fits total.
#--------------------------------------------------------------------------------------
results_all <- list()

for (proportion_sparsity in sparsity_levels) {
  
  sparsity_label <- paste0(proportion_sparsity * 100, "pct_sparsity")
  message("== Simulating data for sparsity condition: ", sparsity_label, " ==")
  
  sim_data <- simulate_condition(proportion_sparsity)
  
  S_t0_i        <- do.call(rbind, sim_data$S_t0_list_all[[1]])
  y0_bacteria_i <- do.call(rbind, sim_data$y0_bacteria_list_all[[1]])
  y_obs_i       <- sim_data$y_obs_all[1, , , ]
  
  data_list_Stan <- list(
    N_met = ncol(S_t0_i),
    N_bac = ncol(y0_bacteria_i),
    T = dim(y_obs_i)[2],
    R = dim(y_obs_i)[1],
    ts = ts_sim,
    y_obs = y_obs_i,
    S_t0 = S_t0_i,
    y0_bacteria = y0_bacteria_i,
    t0 = 0,
    k = k,
    N_produced = length(produced_metabolites),
    N_consumed = length(consumed_metabolites),
    produced_metabolites = produced_metabolites,
    consumed_metabolites = consumed_metabolites,
    C0 = C0,
    use_ar1 = 0,
    delta = delta
  )
  
  # Save the shared data/truth once per sparsity level (identical for all 3 priors)
  saveRDS(sim_data$Cij_list,             paste0("truth_Cij_",  sparsity_label, ".rds"))
  saveRDS(sim_data$y_obs_all,            paste0("data_y_obs_", sparsity_label, ".rds"))
  saveRDS(sim_data$S_t0_list_all,        paste0("S0_",         sparsity_label, ".rds"))
  saveRDS(sim_data$y0_bacteria_list_all, paste0("B0_",         sparsity_label, ".rds"))
  
  for (model_name in prior_names) {
    
    condition_label <- paste0(model_name, "_", sparsity_label)
    message("   -> Fitting: ", condition_label)
    
    fit <- sampling(Stan_models[[model_name]],
                    data = data_list_Stan,
                    iter = 1000,
                    chains = 4,
                    cores = 4,
                    warmup = 500,
                    control = list(adapt_delta = 0.9))
    
    Cij_post       <- rstan::extract(fit, pars = "Cij")$Cij
    Cij_est_median <- apply(Cij_post, c(2, 3), median)
    
    y_sim_post   <- rstan::extract(fit, pars = "y_sim")$y_sim
    log_lik_post <- rstan::extract(fit, pars = "log_lik")$log_lik
    
    saveRDS(fit,            paste0("fit_",           condition_label, ".rds"))
    saveRDS(Cij_est_median, paste0("posterior_Cij_", condition_label, ".rds"))
    saveRDS(y_sim_post,     paste0("y_sim_",          condition_label, ".rds"))
    saveRDS(log_lik_post,   paste0("log_lik_",        condition_label, ".rds"))
    
    results_all[[condition_label]] <- list(
      prior         = model_name,
      sparsity      = proportion_sparsity,
      fit           = fit,
      Cij_true      = sim_data$Cij_list[[1]],
      Cij_est       = Cij_est_median
    )
  }
}

#--------------------------------------------------------------------------------------
# 8. Quick sanity check / comparison across priors within each sparsity level
#    (mean absolute error between true and estimated Cij, per prior x sparsity)
#--------------------------------------------------------------------------------------
comparison_table <- do.call(rbind, lapply(names(results_all), function(nm) {
  res <- results_all[[nm]]
  data.frame(
    condition = nm,
    prior     = res$prior,
    sparsity  = res$sparsity,
    mae       = mean(abs(res$Cij_true - res$Cij_est))
  )
}))

comparison_table <- comparison_table[order(comparison_table$sparsity, comparison_table$prior), ]
print(comparison_table)
saveRDS(comparison_table, "comparison_table_all_conditions.rds")