
#load packages needed
library(deSolve) #solves the ordinary differential equations
library(rstan)   #run Bayesian model written in stan
#data manipulation, visualisation, plots
library(tidyr)
library(ggplot2)
library(dplyr)
library(posterior)
library(stringr)

# Prepare for running stan optimally
#configure R/stan to use available CPU resources to speed up Bayesian sampling
options(mc.cores = parallel::detectCores()) #detect how many CPUs you have to use
rstan_options(threads_per_chain = 2)  #each mcmc uses 2 chains

#--------------------------------------------------------------------------------------
# Generating a bacterial community including metabolic interaction and dynamics 
# in silico using a set of coupled ordinary differential equations
#--------------------------------------------------------------------------------------
# Simulate metabolic and bacterial dynamics using a set of ODEs 


  bacteria_metabolite_ode <- function(t, state, parameters)#define how bacteria and metabolite change over time (rate of change)
    {
    with(as.list(c(state, parameters)), {  #extracts parameters, easier to refer to by name
      #determines bacteria and metabolite number
      N_met <- length(S_t0)
      N_bac <- length(y0_bacteria)
      B_total <- sum(state[(N_met + 1):(N_met + N_bac)]) #sums up whats in those positions
      #eg if metabolite 10 and bacteria 3 then 13 positions but doesnt mean adds up to 13
      #total bacteria abundance
      
      #this creates a vector to store the rate of change of each bacterium and metabolites
      dydt <- numeric(N_met + N_bac)
      
      # Metabolite dynamics (produced or consumed by bacteria)
      for (i in 1:N_met) {  #go through each metabolite one at a time 
        dydt[i] <- 0  #resets the effect for each metabolite
        for (j in 1:N_bac) {  #for each metabolite we go through each bacteria
          dydt[i] <- dydt[i] - Cij[i, j] * state[i] * state[N_met + j]
        }
        #state[i] is the abundance of metabolite i, state[N_met + j] is asking for abundance of bacterium j
        #how does bacterium j affect metabolite i
        #formula says interaction strength * metabolite abundance * bacterium abundance
        #-ci,j means the bacterium decrease the metabolites via consumption
        #FOR EACH METABOLITE - SUBTRACT THE EFFECT OF EACH BACTERIUM 
      }
      
      # Bacteria dynamics (Only consumed metabolites contribute to bacterial growth)
      for (j in 1:N_bac) {   #loops through each bacteria
        growth <- -d_bac  #before anything is added, bacteria has baseline death rate(dying)
        for (i in 1:N_met) { #now every bacteria loops through every metabolite
          growth <- growth + C0 * Cij[i, j] * state[i]  #C0 is scaling constant, cij is interaction
        }
        dydt[N_met + j] <- state[N_met + j] * growth * (k - B_total) #this only happens after every metabolite for the single bacteria 
      } #(k - B_total)= looks at carrying capacity. small total bacteria- bacteria can grow. large bacteria approach k and growth slows down.
        #
      list(dydt) #return the rate of change for the ODE
    })
  }
  
  # Constants
  N_met <-10  # Number of metabolites
  N_bac <- 3   # Number of bacteria
  Reps <- 6     # Number of "experimental" replicates
  Obs <- 4      # Number of "experimental" observations #observations PER REP/time point
  k <- 14       # Carrying capacity
  delta <-1e-6  # Small number to ensure numerical stability
  C0 <- 0.1     # Metabolism conversion coefficient
  d_bac <- 0.1  # Death rate per bacterium
  n_sims <- 4  # Number of independent simulations, each simulation uses a separate set of Cij values
  t0 <- 0       # Starting point time series
  tmax <- 4     # End point time series
  ts_sim <- seq(1e-6, tmax, length.out = Obs)
  #time points ODE will be evaluated
  # Time steps for the ODEs, creates exactly Obs equally spaced time points
  set.seed(42)  #makes random simulations reproducible
  
  # Noise addition
  add_noise <- FALSE      # If TRUE adds random noise
  #true underlying abundance is different to models simulated abundance = noise (even if the numbers are close)
  
             
  Sigma_met <- runif(N_met, min = 0.1, max = 0.5)   #generate random noise size for EACH METABOLITE between 0.1 and 0.5
  Sigma_bac <- runif(N_bac, min = 0.01, max = 0.05) #generate random noise size for EACH BACTERIA between 0.01 and 0.05
  Sigma_noise <- c(Sigma_met, Sigma_bac) #like concantenates noise them [1,2 ,3]+ [4,5 , 6] = [1,2,3,4,5,6]
  #the reason we look at noise is because it makes the simulated data look like real experiments
  #to test model we keep noise off or keep noise to see model figure out interaction with noise
  
  
  
  
  
  # Determine the amount of metabolites per category, see N_met_sim
  # Make sure that the order of metabolites is consumed --> produced
  consumed_metabolites <- c(1:7) #positive
  produced_metabolites <- c(8:10) #negative
  
  # Initialize storage objects
  #n_sims = HOW MANY SEPARATE SIMULATIONS ARE WE GOING TO GENERATE
  #EACH SIMULATION WILL COMPUTE EACH OF THESE
  
  Cij_list <- vector("list", n_sims) #stores true cij matrix for every simulation
  S_t0_list_all <- vector("list", n_sims) #store initial metabolite concentrations
  y0_bacteria_list_all <- vector("list", n_sims) #store  initial bacterial 
  y_obs_all <- array(0, dim = c(n_sims, Reps, Obs, N_met + N_bac))#creates 4 dimensional array containing all simulated observations
  
  for (sim in 1:n_sims) { #run through each simulation to generate true cij values
    
    Cij_consumed <- matrix(runif(length(consumed_metabolites) * N_bac, 0, 1),
    #length(consumed_meta) = number, n_bac = number, runif generates random cij values between 0 and 1, matrix arranges into matrix                       
    #generate positive cij values for metabolites consumed
    #generate negative cij values for metabolites produced
    
                           nrow = length(consumed_metabolites), ncol = N_bac)#tells you size of matrix
    Cij_produced <- matrix(runif(length(produced_metabolites) * N_bac, -1, 0),
    #length(produced meta) = number, n_bac = number, runif generates -cij values between -1 and 0, matrix creates matrix
                           nrow = length(produced_metabolites), ncol = N_bac)#tells size of matrix
    Cij <- matrix(0, N_met, N_bac) #create an interaction matrix with initially zeros
    
    
    Cij[consumed_metabolites, ] <- Cij_consumed 
    #for all the columns(,) take the rows with the consumed metabolites(1-7) and fill them with the positive Cij in the matrix.
    Cij[produced_metabolites, ] <- Cij_produced
    ##for all the columns(,) take the rows with the produced metabolites(8-10) and fill them with the negative Cij in the matrix.
    
    # Save matrix current simulation
    Cij_list[[sim]] <- Cij
    
    #generate the starting metabolite for each replicate
    S_t0_list <- replicate(Reps, {
      #c combines these 2 into 1 vector
      c(runif(length(consumed_metabolites), min = 0.5, max = 1),  #generating 7 consumed metabolites each between 0.5 and 1  
        runif(length(produced_metabolites), min = 1e-3, max = 1e-2)) #generating 3 produced metabolites each between 0.001 and 0.01
    }, simplify = FALSE) #keep each replicate as a separate list
    
    #generate the starting bacteria for each replicate
    y0_bacteria_list <- replicate(Reps,
                                  runif(N_bac, min = 0.05, max = 0.2), #generating random N_bac between 0.05 and 0.2 
                                  simplify = FALSE)#also keep each replication as a separate list
    
    # Save initial conditions
    S_t0_list_all[[sim]] <- S_t0_list #sim will tell you what simulation your on
    y0_bacteria_list_all[[sim]] <- y0_bacteria_list
    
    # Solve the set of ODEs 
    for (r in 1:Reps) {  #go through each replicate one at a time
      #holds starting values
      state_init <- c(S_t0_list[[r]], y0_bacteria_list[[r]]) #initial state is created by combing intial states of bac and meta
      
      parameters <- list( #holds everything ODE needs and tells it how to behave. done via deSolve
        Cij = Cij, #interaction strengths
        k = k, #carrying capacity
        S_t0 = S_t0_list[[r]], #gives the initial meta conc FOR CURRENT REPLICATION
        y0_bacteria = y0_bacteria_list[[r]], #gives initial bac conc for CURRENT REPLICATION
        C0 = C0, #metabolite conversion coefficient
        consumed_metabolites = consumed_metabolites, #tells which metabolites are consumed (1-7)
        produced_metabolites = produced_metabolites #tells which metabolites are produced (8-10)
      )
      
      
      #where bacterial and metabolite dynamics are simulated
      #given the initial conditions and interaction coefficients, simulate how the system changes over time
      out <- ode(y = state_init, times = ts_sim, func = bacteria_metabolite_ode,
                 parms = parameters, method = "ode45")
      
      for (i in 1:(N_met + N_bac)) { #loops through every bacteria and metabolite
        true_vals <- out[, i + 1] #out contains results of ODE solver
        #take all the timepoints from coloumn i+1
        #use +1 to get bacteria+metabolite data, first column is time
        #true_vals is the true predicted abundance of metabolite across all time points
        
        
        if (!add_noise) { #not add noise
          y_obs_all[sim, r, , i] <- true_vals + delta #OBSERVED DATA= ODE PREDICTION+ TINY ADJUSTMENT
          
        } else{
          noise <- rnorm(1, mean = 0, sd = Sigma_noise[i])
        
      
          y_obs_all[sim, r, , i] <- exp(log(true_vals + delta) + noise)
          #combines true ODE values with noise and stores the resulting observations
          
          #SO WHAT COMES OUT OF SIMULATIONS:
          #TRUE INTERACTION MATRICES: Cij_list
          #INITIAL METABOLITE CONCENTRATIONS: S_t0_list_all
          #INITIAL BACTERIAL ABUNDANCES: y0_bacteria_list_all
          #SIMULATED OBSERVATIONS:  y_obs_all
        }
      }
    }
  }
    
  
  # Output:
  # Cij_list:             list of n_sims Cij matrices
  # S_t0_list_all:        list of n_sims lists, each containing 6 metabolite vectors corresponding to the intial values
  # y0_bacteria_list_all: list of n_sims lists, each containing 6 bacterial vectors corresponsing to the intial values
  # y_obs_all:            4D array [n_sims simulations, 6 replicates, Obs, N_met + N_bac]
  
  # Make sure to save the generated observations and variables before parameter 
  # estimation using the Bayesian inference approach
  
  #save(Cij_list, S_t0_list_all, y0_bacteria_list_all, y_obs_all, file = "...")}

#--------------------------------------------------------------------------------------
#Bayesian inference framework
  
#move from simulation to parameter estimation
#--------------------------------------------------------------------------------------
# The Bayesian inference framework is written in Rstan
 Stan_code <- "   #stan code
functions {
  // ODE system describing bacteria-metabolite interactions
  //recreating bac/meta equations in rstan so stan can solve ODE while estimating parameters
  vector Bac_Met(real t, //current time
               vector y, //bac+meta abundancies
                matrix Cij, //interaction strengths
               real k, //carrying capacity 
               vector S_t0, //
               int[] consumed_metabolites, //tells model which metabolites are consumed
               real C0) { //parameter controlling the effect of metabolite consumption on bacterial growth.
    int N_met = rows(Cij); //number of metabolites from rows                               
    int N_bac = cols(Cij);   //number of bacteria from columns                               
    vector[N_met + N_bac] dydt;      //stores rate of change of each variable                      
    real B_total = sum(y[(N_met + 1):(N_met + N_bac)]); //calculates total bacterial abundance in same way   
   
    // Metabolite dynamics
    for (i in 1:N_met) { //go through every metabolite one at a time
      dydt[i] = 0; //initial rate of change 0
        for (j in 1:N_bac) { //for every meta go through every bacteria
          dydt[i] -= Cij[i, j] * y[i] * y[N_met + j];// calculates how bacterium j affects metabolite i
          // interaction strenght * metabolite abundance * bacterial abundance
          //-= means subtract this effect from the dydt[i] which means bacteria consume metabolites(decrease)
      }
    }

    // Bacterial growth dynamics
    for (j in 1:N_bac) { //go through bacteria 1 at a time
      real growth = -0.1; //baseline loss/decline rate for the bacteria before considering food from consumed metabolites.
        for (i in consumed_metabolites) { //NOW ONLY LOOK AT METABOLITES THAT BACTERIA CONSUME
          growth += C0 * Cij[i, j] * y[i]; //Add the contribution of each consumed metabolite to bacterial growth.
          //the += take exiting growth and add this contribution to it 
  }
  dydt[N_met + j] = y[N_met + j] * growth * (k - B_total); //rate of change of bacterium j
  //bacterial abundance * growth * how much room left before reaching carrying capacity
}
    return dydt; //returns rates of change to ODE solver
  }
}

data { //defines everything stan receives as known information
  int<lower=1> R; //replicates
  int<lower=1> T; //timepoints
  int<lower=1> N_met; // metabolites
  int<lower=1> N_bac;//bacteria
  array[T] real ts; //T = number of timepoints, real=can have decimals, ts = timepoints when observations r taken
  array[R, T] vector[N_met + N_bac] y_obs; //y-obs[2,5,12] = 2 reps, time = 5, bacterium 2
  vector[N_met] S_t0[R]; //stores the initial metabolite conc for each rep
  vector[N_bac] y0_bacteria[R]; //stores inital bacteria conc for each rep
  real t0; //starting time of ODE
  real<lower=0> k; //k must be bigger than 0 to be positive 
  int<lower=0> N_produced; //cant be negative
  int<lower=0> N_consumed; //cant be negative
  int produced_metabolites[N_produced]; //stores indices of produced metabolites 
  int consumed_metabolites[N_consumed];  //stores indices of consumed metabolites
  int<lower=0, upper=1> use_ar1; //on/of switch for AR1. //0= no noise; 1=noise
  real<lower=0> delta; //delta cant be 0
  real<lower=0, upper=1> C0; //effect of consumed metabolites on bacteria
}

parameters {
  // Cij matrix, the priors for the cij matrix are defined here since those are uniform and not in the model block
  // Make sure the order of the metabolites is consumed --> produced
  
  matrix<lower=0, upper=1>[N_consumed, N_bac] Cij_pos;  // Uniform prior between 0 and 1
  //cij pos=represents positive interactions for consumed metabolites between 0 and 1
  //[row,column]
  
  matrix<lower=-1, upper=0>[N_produced, N_bac] Cij_neg; // Uniform prior between -1 and 0
  //cij neg = represents negative interactions for produced metabolites between -1 and 0
  //[row,column]
  real<lower=0> sigma_ar; }
  
transformed parameters{
   

  // Construct full Cij matrix by combining positive and negative parts
  matrix[N_met, N_bac] Cij; //creates empty cij matrix: row meta columns bac
  
  for (i in 1:N_consumed) { //loop through consumed meta
    for (j in 1:N_bac) { //loop thro bac
      Cij[i, j] = Cij_pos[i, j]; //copies the positive interaction value into full cij matrix
    }
  }
  
  for (i in 1:N_produced) {  //loop through produced meta
    int idx = N_consumed + i; //this is to find correct row where produced meta starts
    for (j in 1:N_bac) { //loop through bac
      Cij[idx, j] = Cij_neg[i, j]; //copies the negative interaction value into full matrix
    }
  }
}

//simulate system over time for each replicate
model { //model calculations begin

  array[R, T] vector[N_met + N_bac] mu; //mu created to store the ODE PREDICTED values 
  //R replicates; T time points; meta + bac conc
  //mu[3,4,13] = 3 replicates, time =5 and bacterium 3
  vector[N_met + N_bac] y0; //create vector to hold inital state of system




  for (r in 1:R) {  //loop through each replicate 
    y0[1:N_met] = S_t0[r];  //Puts the initial metabolite concentrations for replicate r into the first part of y0.
    y0[(N_met + 1):(N_met + N_bac)] = y0_bacteria[r]; //Puts the initial bacterial abundances into the rest of y0.
    
    //run the Bac_Met ODE using these starting values and parameters, and give me the predicted values at all the specified time points.
    mu[r] = ode_rk45(Bac_Met, //ODE function that has equations for bac and meta
                     y0,//initial concentrations
                     t0,//initial time
                     ts, //time we want predictions
                     Cij,//interaction strength
                     k,
                     S_t0[r], //initial metabolites
                     consumed_metabolites,
                     C0);
  }






// Prior for the residual noise
sigma_ar ~ normal(0, 0.2); //Gives the noise size a prior centred around 0 with SD 0.2, how much random variation there is around the ODE predictions.

 
for (r in 1:R) { //loop through every replicate 
    for (l in 1:(N_met + N_bac)) { //go through every meta and bac
      vector[T] eps; //create vector to store error at each time point
      for (t in 1:T) { //go through each time point
        eps[t] = log(y_obs[r, t, l] + delta) - log(mu[r, t, l] + delta);
        //in log scale calculate observed value- model predicted value
        //if model observation very close to prediction esp = 0
      }
      eps ~ normal(0, sigma_ar);
      
      
         
        }
      }
    }
  

//dont run yet
//generated quantities { //used for prediction/diagnostics
  //array[R, T] vector[N_met + N_bac] y_sim; //y_sim stores the ODE predicted state, R=replicates, T=time points, N_met+N_bac = total variables
  
  //vector[N_met + N_bac] y0_sim; //stores initial conditions of ODE
  //array[R, T, N_met + N_bac] real log_lik; //creates 3D array to store log-likelihood contribution for every r,t,species
  //real log_lik is the contribution from that particular observation
  
  
  //real rho_pop = 2 * inv_logit(mu_rho) - 1;   //creates transformed version of mu_rho
  //inv_logit converts mu_rho into interval 0,1 then 2*()-1 = maps it (-1,1)
  //this is valid range for AR1 correlation paramter represented by rho_pop

  //for (r in 1:R) {  //loops over replicate
    //y0_sim[1:N_met] = S_t0[r]; //containes initial metabolite state for each replication r
    //y0_sim[(N_met + 1):(N_met + N_bac)] = y0_bacteria[r]; //fill the rest of y0_sim with bacterial abundances
    
    //colve ODE system 
    //y_sim[r] = ode_rk45(Bac_Met,    
                        //y0_sim, //provides solver initial state
                        //t0,
                        //ts,
                        //Cij,
                        //k,
                        //S_t0[r], //initial metabolite conc
                        //consumed_metabolites,
                        //C0);


    //for (l in 1:(N_met + N_bac)) { //loop through every state
      //for (t in 1:T) { //loop through every observation of time
        //real eps = log(y_obs[r, t, l] + delta) - log(y_sim[r, t, l] + delta); //calculates the log-scale residual between observation and model prediction
        //if (use_ar1 == 1) { //residuals correlated
          //if (t == 1) { //checks if this is first time point
          
            //log_lik[r, t, l] = normal_lpdf(eps | 0, sigma_ar / sqrt(1 - square(rho[r]))); //calculates the log probability of the first residual
         // } else { //for every time point after the first
          
            //real prev_eps = log(y_obs[r, t - 1, l] + delta) - log(y_sim[r, t - 1, l] + delta); //calculate residual from previous timepoint
            //log_lik[r, t, l] = normal_lpdf(eps | rho[r] * prev_eps, sigma_ar); //actual Ar1 relationship
          //}
       // } else { //when model not using AR1
         // log_lik[r, t, l] = normal_lpdf(eps | 0, sigma_ar); //independent normally distributed residuals, no dependence on previous time point
       // }
      //}
    //}
  //}
//}
"

# Compiling and parameter estimation is written in R
# Make sure to have saved the input variable data, the datalist is not saved here
# Compile the model

#takes your Stan code stored as text in Stan_code and compiles it into a Stan model object that can be passed to sampling()
Stan_code <- stan_model(model_code = Stan_code)

 # Storage objects
fits_test <- vector("list", n_sims)  #create an empty list with n_sims slots to store the Stan fit objects.
true_Cij_list <- vector("list", n_sims) #create an empty list to store the true Cij matrices used to generate each simulated dataset.
posterior_Cij_list <- vector("list", n_sims) #creates an empty list to store the estimated/posterior Cij matrices from Stan
n_sims <- 4 # Number of simulations

# Parameter estimation
for (i in 1:n_sims) { #loop over simulations
  S_t0_i <- do.call(rbind, S_t0_list_all[[i]]) #initial metabolite conc for simulation i, combines elements of list using rbind
  y0_bacteria_i <- do.call(rbind, y0_bacteria_list_all[[i]]) #contains initial bacterial abundances for each rep
  y_obs_i <- y_obs_all[i, , , ] #extracts the observed data corresponding to simulation i
  #select simulation i, select t , select species, select r
  
  # Create a data list for all data Stan expects
  data_list_Stan <- list(
    N_met = ncol(S_t0_i), #counts number of metabolite columns
    N_bac = ncol(y0_bacteria_i), #counts number of bacterial columns
    T = dim(y_obs_i)[2], #gets size of second dimension- t
    R = dim(y_obs_i)[1], #gets size of replicates
    ts = ts_sim, #observation time points
    y_obs = y_obs_i, #pass the observed metabolite/bacteria measurements to Stan
    S_t0 = S_t0_i, #pass the initial metabolite concentrations to Stan.
    y0_bacteria = y0_bacteria_i, #pass the initial bacterial abundances to Stan
    t0 = 0, #initial ODE time to zero
    k = 14,
    N_produced = length(produced_metabolites), #count metabolites produced
    N_consumed = length(consumed_metabolites), #count metabolites consumed
    produced_metabolites = produced_metabolites, #pass the list identifying the produced metabolites to Stan
    consumed_metabolites = consumed_metabolites, #pass the list identifying the consumed metabolites to Stan
    C0 = 0.1,
    use_ar1 = 0, #tells Stan not to use the AR1 residual structure
    delta = delta #passes the R variable delta into Stan
  )
  
  # Fit Stan model 
  fit <- sampling(Stan_code, #run Bayesian sampling using the compiled Stan mode -stored in fit
                  data = data_list_Stan, #provides Stan with data  just constructed
                  iter =1000,    #4000
                  chains =1,    #4  #independent mcmc chains
                  cores = 2,    #4 CPU cores to be used for the sampling
                  warmup =500, #2000  #warmup iterations used to tune sampler and not used for posterior inference
                  control = list(adapt_delta = 0.95)) #this makes sure not too conservative
  
  # Save the fit object
  fits_test[[i]] <- fit  #stores fitted Stan model for simulation i in list
  
  # Save true Cij values used to generate this simulation i
  true_Cij_list[[i]] <- Cij_list[[i]]  
  
  # Extract posterior median for Cij
  #extracts all posterior samples of the parameter Cij from the fitted Stan model
  #So Cij_post contains the posterior distribution of every element of Cij
  Cij_post <- rstan::extract(fit, pars = "Cij")$Cij
  posterior_Cij_list[[i]] <- apply(Cij_post, c(2,3), median) #take median across MCMC samples for every element of the Cij matrix
  #keep dimensions 2 and 3 and calculate median over remaining dimension
}

 

saveRDS(fit, file = ".......rds")   #saves fit
saveRDS(true_Cij_list, file = "....rds") #saves  true cij
saveRDS(posterior_Cij_list, file = "....rds") #saves posterior median cij estimates
