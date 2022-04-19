# DESCRIPTION -------------------------------------------------------------

# This program estimates the RDEU model with CRRA and hyperbolic discounting based on Multiple Price Lists decision problems.
# Reference data: Andersen et al. (2008)



# PATHS AND PACKAGES ------------------------------------------------------

# Clean up environment
rm(list=ls())

# Load required packages
packages_required <- c("MASS", "pracma", "msm", "here")
install.packages( setdiff( packages_required, rownames( installed.packages() ) ) )
for(name in packages_required) { library(name, character.only = TRUE) }

# Set paths
path                = here::here()
path = paste0(path, "/Desktop/random-models-joint-risk-time")
path_empirical_data = paste0(path, "/Data/processed_data/AHLR/")



# GENERAL SETTINGS --------------------------------------------------------

# SET: Level of background consumption
# Andersen et al. (2008) use omega = 118 
# Alternative specification with omega = 3304 (not reported in paper)
omega = 118

# Utility specification
util_spec = "CRRA"



# EXPERIMENTAL PARAMETERS -------------------------------------------------

# RISK TASKS 

# Payments associated with the safe and risky lottery
safe_lotteries  = list(c(2000, 1600), 
                       c(2250, 1500),
                       c(2000, 1750),
                       c(2500, 1000))
risky_lotteries = list(c(3850,  100),
                       c(4000,  500),
                       c(4000,  150),
                       c(4500,   50))
# Probability sequence (probabilities of being in the good state)
# NOTE: We omit p=1 here, since it has no relevance for risk aversion
p_seq = seq(0.1, 0.9, 0.1)


# TIME TASKS

# Payment dates
t_0 = 1                       # Earlier payment date (in months)
t_1 = c(2, 5, 7, 13, 19, 25)  # Later payment dates  (in months) [Andersen et al. (2008), p.588]

# Payment amounts
x_0 = rep(3000, 10)           # Earlier payment list 
later_payments_list = list(   # Later payment lists
  "x_2" =c(3012, 3025, 3037, 3049, 3061, 3073, 3085, 3097, 3109, 3120),
  "x_5" =c(3050, 3100, 3151, 3202, 3253, 3304, 3355, 3407, 3458, 3510),
  "x_7" =c(3075, 3152, 3229, 3308, 3387, 3467, 3548, 3630, 3713, 3797),
  "x_13"=c(3153, 3311, 3476, 3647, 3823, 4006, 4196, 4392, 4595, 4805),
  "x_19"=c(3232, 3479, 3742, 4020, 4316, 4630, 4962, 5315, 5687, 6082),
  "x_25"=c(3313, 3655, 4027, 4432, 4873, 5350, 5869, 6431, 7039, 7697)
)

# Initialize grid of r-values
# NOTE: r_grid needs to be dense "enough" in the sense that estimated parameters
# are not sensitive to a change in the r_grid.
r_grid = seq(-30, 30, 1/500)



# AUXILIARY FUNCTIONS -----------------------------------------------------

# DEF FUNCTION: Calculate conditional (on r) location parameter
# Return: Conditional location parameter
mu_gamma_cond_r = function(r, mu_r, mu_gamma, sigma_r, sigma_gamma, rho) {
  
  mu_cond = mu_gamma + (sigma_r/sigma_gamma) * rho * (r - mu_r)
  
  return (mu_cond)
}


# DEF FUNCTION: Calculate conditional (on r) scale parameter
# Return: Conditional scale parameter
sigma_gamma_cond_r = function(r, sigma_gamma, rho){
  
  sigma_cond = sqrt(1-rho^2) * sigma_gamma
  
  return (rep(sigma_cond, length(r)))
}



# SPECIFY: UTILITY FUNCTIONS ----------------------------------------------

# DEFINE: CRRA utility function
utility_function = function(r, omega, x){
  
  if (r != 1){
    util_ = { (omega + x)^(1-r) - omega^(1-r) }/{ 1-r }
  }
  if (r == 1){
    util_ = log(omega + x) - log(omega)
  }
  
  return (util_)
}



# THRESHOLD FUNCTIONS -----------------------------------------------------

# DEF FUNCTION: Auxiliary function to calculate risk thresholds
# Return: Auxiliary function to facilitate calculation of risk thresholds
create_aux_r_threshold_fun = function(r, p, omega, x_s, x_r){
  
  aux_r_threshold = function(r){
    # Define the threshold (as given by the formula in the paper)
    r_threshold = p * ( utility_function(r, omega, x_s[1]) - utility_function(r, omega, x_r[1]) ) +
      (1-p) * ( utility_function(r, omega, x_s[2]) - utility_function(r, omega, x_r[2]) )
    
    return (r_threshold)
  } 
  
  return (aux_r_threshold)
}



# DEF FUNCTION: Calculate gamma thresholds, conditional on r-value
# Return: Gamma threshold
gamma_threshold_fun = function(r, omega, x_0, x_1, t_0, t_1){
  
  u_0 = utility_function(r, omega, x_0)
  u_1 = utility_function(r, omega, x_1)
  
  if (t_1*u_0 <= t_0*u_1){
    d = 0
  }  
  else if (t_1*u_0 > t_0*u_1){
    d = (u_1 - u_0)/(t_1*u_0 - t_0*u_1)
  } 
  
  g = log(d)
               
  return (g)
}

# NOTE: Threshold function is different due to hyperbolic discounting
# See Appendix F. for the threshold



# CALCULATION OF THRESHOLDS -----------------------------------------------

# DEF FUNCTION: Create r threshold matrix
# Return: Matrix that contains risk thresholds
#         rows - value of p of p_seq
#         cols - payment combination
create_r_threshold_matrix = function(omega){
  
  # Initialize empty matrix
  threshold_matrix_risk = matrix(
    data = NA,
    nrow = length(p_seq),
    ncol = length(risky_lotteries)
  )
  
  for (p_ in 1:length(p_seq)){  # Iterate over the probabilities
    for (i_ in 1:length(risky_lotteries)){  # Iterate over the payment combinations
      
      aux_r_threshold_func = 
        # Create the auxiliary threshold function (a function of r only)
        create_aux_r_threshold_fun(
          r, p_seq[p_], omega, 
          safe_lotteries[[i_]], 
          risky_lotteries[[i_]]
        )
      
      # Find the unit roots of the threshold function
      value_r_threshold = uniroot(aux_r_threshold_func, c(-20, 20))$root 
      
      # Store threshold for a given (probability, payment) combination
      threshold_matrix_risk[p_, i_] = value_r_threshold
      
    }  # end for i_
  }  # end for p_
  
  return(threshold_matrix_risk)
}

# CREATE: Threshold matrix risk
threshold_matrix_risk = create_r_threshold_matrix(omega)


# DEF FUNCTION: Populate matrix for gamma-thresholds, conditional on a given r-value
# Return: Matrix that contains gamma thresholds, conditional on a given r-value
#         rows - value of payment list
#         cols - timing of payment
create_gamma_threshold_matrix = function(r){
  
  # Initialize empty matrix to store gamma thresholds
  gamma_threshold_matrix = matrix(
    data = NA, 
    nrow = length(later_payments_list[[1]]),  # Each row corresponds to a rate of return
    ncol = length(t_1)                        # Each col corresponds to a time gap
  )
  
  for (t_ in 1:length(t_1)){
    for (x_1_ in 1:length(later_payments_list[[t_]])){
      
      value_ = gamma_threshold_fun(
        r, omega,
        x_0[x_1_], later_payments_list[[t_]][x_1_], 
        t_0, t_1[t_]
      )
      gamma_threshold_matrix[x_1_, t_] = value_
    }
  }
  return (gamma_threshold_matrix)
}

# NOTE: Function adjusted to accommodate hyperbolic discounting case


# CREATE: List of gamma-threshold matrices for every r-value of r_grid
list_gamma_threshold_matrices = list()
for (r_ in 1:length(r_grid)){
  list_gamma_threshold_matrices[[r_]] = create_gamma_threshold_matrix(r_grid[r_])
  # NOTE: Additional condition in hyperbolic discounting to ensure that thresholds != +/- Inf
  list_gamma_threshold_matrices[[r_]][list_gamma_threshold_matrices[[r_]]==-Inf] = -99999999999
}



# CALCULATION OF CHOICE PROBABILITIES -------------------------------------

# DEF FUNCTION: Calculate marginal probabilities
# Return: Grid of marginal probabilities
calc_marginal_probabilities = function(pars){

  # Retrieve the parameters
  if (length(pars) > 2){  
    mu_r = pars[1]
    sigma_r = pars[3]
  }
  if (length(pars) == 2){ 
    mu_r = pars[1]
    sigma_r = pars[2]
  }

  # Calculate the marginal probabilities (by proper normalization)
  value_ = dnorm(r_grid, mu_r, sigma_r) / sum(dnorm(r_grid, mu_r, sigma_r))

  return(value_)
}


# DEF FUNCTION: Calculate choice probabilities for risk tasks
# Return: Matrix of choice probabilities for the risk task
#         rows - choice probability for given p
#         cols - choice probability for given payment combination
calc_choice_prob_mat_risk = function(pars){
  
  # Retrieve parameters
  if (length(pars) > 2){  
    mu_r    = pars[1]
    sigma_r = pars[3]
  }
  if (length(pars) == 2){ 
    mu_r    = pars[1]
    sigma_r = pars[2]
  }
  
  # Calculate the choice probabilities
  choice_probs = 1 - pnorm(threshold_matrix_risk, mu_r, sigma_r)
  
  return(choice_probs)
}


# DEF FUNCTION: Calculate conditional (on r) choice probabilities of earlier [a] and later [b] choices
# Return: List of conditional choice probability matrices
#         [[x]][[ ]]: x="a" - a-choices, 
#                     x="b" - b-choices
#         [[ ]][[x]]: matrix of conditional choice probabilities for a given r-value
calc_cond_probs_time = function(pars){
  
  # Initialize list to store conditional choice probability matrices
  list_cond_choice_prob_matrices_a = list()
  list_cond_choice_prob_matrices_b = list()
  
  # Restore parameters
  mu_r        = pars[1]
  mu_gamma    = pars[2]
  sigma_r     = pars[3]
  sigma_gamma = pars[4]
  rho         = pars[5]
  
  # Calculate conditional location and scale parameters
  mu_gamma_r    = mu_gamma_cond_r(r_grid, mu_r, mu_gamma, sigma_r, sigma_gamma, rho)
  sigma_gamma_r = sigma_gamma_cond_r(r_grid, sigma_gamma, rho)
  
  # Calculate conditional choice probabilities
  # NOTE: The interpretation of the time parameter is inverse to the exponential case. 
  #       That means, alternative 1 is chosen when the parameter is below the conditional threshold
  #       Consequently, conditional choice probability of alternative 0 is now given by "1 - pnorm(threshold)"
  for (r_ in 1:length(list_gamma_threshold_matrices)){
    list_cond_choice_prob_matrices_a[[r_]] = 1 - pnorm(
      list_gamma_threshold_matrices[[r_]], 
      mu_gamma_r[r_], sigma_gamma_r[r_]
      )
    list_cond_choice_prob_matrices_b[[r_]] = pnorm(
      list_gamma_threshold_matrices[[r_]], 
      mu_gamma_r[r_], sigma_gamma_r[r_]
    )
  }  # end for r_
    
  # Summarize conditional (on r) choice probabilities in a list
  list_cond_choice_prob_matrices = list(
    "a"=list_cond_choice_prob_matrices_a, "b"=list_cond_choice_prob_matrices_b
    )
  
  return(list_cond_choice_prob_matrices)
}


# DEF FUNCTION: Calculate weighted choice probabilities of earlier [a] and later [b] choices
# Return: List of matrices that contain the weighted sum of (conditional on r) choice probabilities ...
#         [[x]][[ ]]: x="a": ... for alternative a
#         [[x]][[ ]]: x="b": ... for alternative b
calc_weighted_sum_cond_probs_time = function(pars){
 
  # Retrieve the conditional choice probabilities for a (earlier) and b (later)
  list_cond_choice_prob_matrices_a = calc_cond_probs_time(pars)[["a"]]
  list_cond_choice_prob_matrices_b = calc_cond_probs_time(pars)[["b"]]
  
  # Calculate the pdf values at every value of r_grid 
  pdf_vals = dnorm(r_grid, pars[1], pars[3])
  
  # Initialize empty matrix to capture the weighted choice probabilities
  weighted_sum_prob_a = matrix(
      0,
      nrow = nrow(list_gamma_threshold_matrices[[1]]),
      ncol = ncol(list_gamma_threshold_matrices[[1]])
    )
  weighted_sum_prob_b = weighted_sum_prob_a
  
  # Sum the weighted choice probabilities (for every value of r_grid)
  for (r_ in 1:length(list_gamma_threshold_matrices)){
    weighted_sum_prob_a = weighted_sum_prob_a + list_cond_choice_prob_matrices_a[[r_]] * pdf_vals[r_]
    weighted_sum_prob_b = weighted_sum_prob_b + list_cond_choice_prob_matrices_b[[r_]] * pdf_vals[r_]
  }  # end for r_
  
  # Store weighted (conditional on r) choice probabilities in a list
  list_weighted_sum_probs = list(
    "a"=weighted_sum_prob_a, "b"=weighted_sum_prob_b
    )

  return(list_weighted_sum_probs)
}


# DEF FUNCTION: Calculate sum of pdf values
# Return: Sum of pdf values over r_grid
calc_sum_pdf_vals = function(pars){
  
  # Retrieve parameters
  if (length(pars) > 2){ 
    mu_r    = pars[1]
    sigma_r = pars[3]
  }
  if (length(pars) == 2){
    mu_r    = pars[1]
    sigma_r = pars[2]
  }
  
  sum_pdf_vals = sum(dnorm(r_grid, mu_r, sigma_r))
  
  return(sum_pdf_vals)
}
  


# DATA: LOAD EMPIRICAL DATA -----------------------------------------------

# Retrieve meta-information (individuals)
num_individuals = length(readRDS(file = paste0(path_empirical_data, "id_list.rds")))
ids_individuals = readRDS(file = paste0(path_empirical_data, "id_list.rds"))

# DEF FUNCTION: Extract empirical data
# Return: List of empirical choice matrices for risk and time
get_empirical_data = function(path_data_risk, path_data_time){
  
  # Read empirical risk data
  readfile_risk                  = readRDS(paste0(path_empirical_data, path_data_risk, ".rds"))
  empirical_matrix_risk_choice_a = readfile_risk[["a"]][1:9, ]  # Discard p=1 
  empirical_matrix_risk_choice_b = readfile_risk[["b"]][1:9, ]  # Discard p=1 
  empirical_choice_matrix_risk   = list(empirical_matrix_risk_choice_a, empirical_matrix_risk_choice_b)
  
  # Read empirical time data
  readfile_time                  = readRDS(paste0(path_empirical_data, path_data_time, ".rds"))
  empirical_matrix_time_choice_a = readfile_time[["a"]]
  empirical_matrix_time_choice_b = readfile_time[["b"]]
  empirical_choice_matrix_time   = list(empirical_matrix_time_choice_a, empirical_matrix_time_choice_b)
  
  empirical_choice_matrices = list(
    empirical_choice_matrix_risk,
    empirical_choice_matrix_time
  )
  
  return(empirical_choice_matrices)
}

# Get empirical data for representative agent case
emp_data = get_empirical_data(
  "AHLR_empirical_matrices_risk", 
  "AHLR_empirical_matrices_time"
  )

# Retrieve empirical data for each individual
list_emp_data_individual = list()
for (id_ in seq_along(ids_individuals)){
  
  list_emp_data_individual[[id_]] = get_empirical_data(
    paste0("individual/AHLR_empirical_matrices_risk_ind_", ids_individuals[id_]), 
    paste0("individual/AHLR_empirical_matrices_time_ind_", ids_individuals[id_])
  )
}


# ESTIMATION: DEFINE LOG-LIKELIHOOD FUNCTION ------------------------------
 
# DEF FUNCTION: Log-likelihood function (risk and time, joint estimation)
# Return: Value of the log-likelihood function
log_likelihood_function = function(pars, data){
  
  # Retrieve choice probabilities for risk
  cp_matrix_risk     = calc_choice_prob_mat_risk(pars)
  log_cp_matrix_risk = log(cp_matrix_risk)
  
  # Assemble log-likelihood function for risk and sum likelihood contributions
  likelihoods_risk = data[[1]][[1]] * log_cp_matrix_risk + data[[1]][[2]] * log(1-cp_matrix_risk)
  likelihood_risk  = sum(likelihoods_risk)
  
  # Set finite value in case the log-likelihood function is not defined
  if (is.na(likelihood_risk) | is.nan(likelihood_risk) ){
    likelihood_risk = -999999999
  }
  
  # Retrieve weighted conditional choice probabilities for time
  log_ws_ccp_a     = log(calc_weighted_sum_cond_probs_time(pars)[["a"]])
  log_ws_ccp_b     = log(calc_weighted_sum_cond_probs_time(pars)[["b"]])
  log_sum_pdf_vals = log(calc_sum_pdf_vals(pars))
  
  # Assemble log-likelihood function for time and sum likelihood contributions
  likelihoods_time = data[[2]][[1]] * (log_ws_ccp_a - log_sum_pdf_vals) + data[[2]][[2]] * (log_ws_ccp_b - log_sum_pdf_vals)
  likelihood_time = sum(likelihoods_time)
  
  # Add up the likelihood contributions from risk and time
  likelihood = likelihood_risk + likelihood_time
  # Set finite value in case the log-likelihood function is not defined
  if (is.nan(likelihood)){
    likelihood = -999999999
  }
  
  print(paste0("Likelihood-value: ", likelihood))
  
  # MLE procedure L-BFGS-B requires finite values in the definition of the log-likelihood function
  if (likelihood == -Inf){
    likelihood = -999999999
  }
  
  return(-likelihood)
}



# ESTIMATION: OPTIMIZATION ------------------------------------------------

# Initialize starting parameters
starting_pars_ = c(0.82, -5.2, 1.06, 1.06, -0.70)

# Set bounds
lower_bounds_  = c(-10, -10, 0.01, 0.01, -0.99) 
upper_bounds_  = c( 10,  10,   10,   10,  0.99)


# Optimization routine
mle = optim(
  par     = starting_pars_, 
  fn      = log_likelihood_function, 
  lower   = lower_bounds_, 
  upper   = upper_bounds_,
  data    = emp_data, 
  method  = "L-BFGS-B", 
  hessian = TRUE
  )



# ESTIMATION: CALCULATION STANDARD ERRORS ---------------------------------
# NOTE: A detailed description is available in the ReadMe

# Extract gradients for each individual at the ML estimates and store them in a list 
list_mult_gradients_individual = list()
for (id_ in seq_along(ids_individuals)){
  
  # DEF FUNCTION: Wrapper for the log-likelihood function 
  # Return: Log-likelihood function as function of the parameters only, given individual choices
  aux_log_likelihood_function = function(pars){
    
    aux_ll_function_ = log_likelihood_function(pars, list_emp_data_individual[[id_]])
    
    return(aux_ll_function_)
  }
  
  # Store the multiplied gradients in a list
  gradient                              = grad(aux_log_likelihood_function, mle$par)
  mult_gradients                        = gradient %*% t(gradient)
  list_mult_gradients_individual[[id_]] = mult_gradients
}

# Sum the gradients and get the inverse of the Hessian
sum_gradients_individual = Reduce("+", list_mult_gradients_individual)
inv_hessian              = solve(mle$hessian)

# Calculate robust standard errors, clustered on individual level
robust_cov_matrix = inv_hessian %*% sum_gradients_individual %*% inv_hessian 
est_se            = sqrt(diag(robust_cov_matrix))


# ESTIMATION: EXPORT PARAMETERS -------------------------------------------

# Retrieve estimated parameters
est_pars = mle$par
# NOTE: Moments of "delta" are not reported in the paper.
# However, "mu_delta" can be obtained by non-linear transformation.
# est_mu_delta = exp(est_pars[2])/(1+exp(est_pars[2]))
# The standard errors of "sigma_delta" can be obtained by a bootstrapping procedure.

# Retrieve grid information
r_grid_info = c(min(r_grid), max(r_grid), round(r_grid[2] - r_grid[1], 5))

# Assemble export vector
export_vector = c(
  "mu_r"             = est_pars[1], "sigma_r"      = est_pars[3], 
  "mu_g"             = est_pars[2], "sigma_g"      = est_pars[4],  
  "rho"              = est_pars[5],
  "se_mu_r"          = est_se[1], "se_sigma_r"     = est_se[3], 
  "se_mu_gamma"      = est_se[2], "se_sigma_gamma" = est_se[4], 
  "se_rho"           = est_se[5], 
  "likelihood-value" = mle$value, "omega" = omega,
  "min_r_grid"       = r_grid_info[1], "max_r_grid" = r_grid_info[2], "step_r_grid" = r_grid_info[3]
)

# Export parameter estimates
write.table(export_vector, 
            file=paste0(path, "/Output/estimates/AHLR/", util_spec, "_", "estimated_pars_risk_and_time_omega_", omega, "_hyperbolic", ".csv"),
            row.names=TRUE, col.names=FALSE, sep=",")
saveRDS(export_vector, file=paste0(path, "/Output/estimates/AHLR/", util_spec, "_", "estimated_pars_risk_and_time_omega_", omega, "_hyperbolic", ".rds"))
