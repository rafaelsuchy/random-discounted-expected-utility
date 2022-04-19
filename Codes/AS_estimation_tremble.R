# DESCRIPTION -------------------------------------------------------------

# This program estimates the RDEU model with tremble based on Convex Menus decision problems.
# Reference data: Andreoni and Sprenger (2012)



# PATH AND PACKAGES ---------------------------------------------------

# Clean up environment
rm(list = ls())

# Load required packages
packages_required <- c("MASS", "pracma", "msm", "here")
install.packages( setdiff( packages_required, rownames( installed.packages() ) ) )
for(name in packages_required) { library(name, character.only = TRUE) }

# Set paths
path                = here::here()
path_empirical_data = paste0(path, "/Data/processed_data/AS/")



# GENERAL SETTINGS --------------------------------------------------------

# SET: Level of background consumption
# Andreoni and Sprenger (2012, Appendix p.6) use 4.05 (as "average")
# Alternative specification with omega = 113.40 (not reported in paper)
omega = 4.05

# Initialize grid of r-values
# NOTE: r_grid needs to be dense "enough" in the sense that estimated parameters
# are not sensitive to a change in the r_grid.
r_grid = seq(-50, 50, 1/500)

# Initialize a-threshold values
a_threshold_values = c(0.04, 0.14, 0.24, 0.34, 0.44, 0.54, 0.64, 0.74, 0.84, 0.94, 1)
# Initialize a-shares
a_shares = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)



# EXPERIMENTAL PARAMETERS -------------------------------------------------

# Probability sequences
p_0 = c(1.0, 1.0, 0.8, 0.5, 0.5, 0.4)
p_1 = c(1.0, 0.8, 1.0, 0.5, 0.4, 0.5)

# Timings
t_0 = c(7)           # after 7 days
t_1 = c(28+7, 56+7)  # after 35 / 63 days 

# Re-scale timing to months
t_0 = t_0 / 28  
t_1 = t_1 / 28

# Payment amounts
payments_7_35 = list(  # payments after 7 and 35 days
  c(20, 19, 18, 17, 16, 15, 14), 
  c(20, 20, 20, 20, 20, 20, 20)
)
payments_7_63 = list(  # payments after 7 and 63 days
  c(20, 19, 18, 17, 16, 15, 14), 
  c(20, 20, 20, 20, 20, 20, 20)
)
list_payments = list(payments_7_35, 
                     payments_7_63)

# Structure of list_payments
#   [[x]][[ ]]: index for (t_0, t_1) combination 
#               x=1: (7, 35) [short time gap]
#               x=2: (7, 63) [long time gap]
#   [[ ]][[y]]: index for payment combinations
#               y=1: t_0 payments
#               y=2: t_1 payments



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


# DEF FUNCTION: Set every input larger 1 to 1
# Return: Input or 1
# Reason: Any threshold delta >= 1 will be set to 1
min_one_value = function(eval_value){
  
  if (eval_value>=1){
    eval_value = 1
  }
  
  return(eval_value)
}


# DEF FUNCTION: Sum pdf-values
# Return: Sum of pdf values for any r-grid point
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



# CALCULATION OF THRESHOLDS -----------------------------------------------

# DEF FUNCTION: Calculate gamma thresholds as function of experimental parameters, given r and a values
# Return: Gamma thresholds
# NOTE: The function reflects the equation in the paper
calc_gamma_threshold = function(p_0, p_1, r, omega, a, x_0, x_1, t_0, t_1){
  
  delta = c()
  
  if (r <= 0){ # Convex monetary utility: unique threshold, independent of a --> only corner solutions have non-null probability
    delta_aux = {p_0/p_1}^{1/(t_1 - t_0)} * { ( (omega + x_0)^(1-r) - omega^(1-r) )/( (omega + x_1)^(1-r) - omega^(1-r) ) }^{1/(t_1 - t_0)}
    delta = sapply(delta_aux, min_one_value)
  }
  
  if (r > 0){ # Concave monetary utility: interior solutions possible
    delta_aux = { (p_0/p_1) * (x_0/x_1) * ( (omega + a * x_1)/(omega + (1-a)*x_0) )^r }^{1/(t_1 - t_0)}
    delta = sapply(delta_aux, min_one_value)
  }
  
  # Obtain gamma threshold through transformation (as described in the paper)
  gamma = log(delta/(1-delta))
  
  return(gamma)
}



# DEF FUNCTION: Initialize gamma threshold matrix
# Return: Matrix with #payments rows and #r_grid-points columns
init_matrix_r_grid = function(list_payments, r_grid){
  
  data = matrix(
    data = NA, 
    nrow = length(list_payments[[1]][[1]]),
    ncol = length(r_grid)
  )
  
}


# DEF FUNCTION: Create threshold matrices
# Structure indices:
#                   [[x]][[ ]][[ ]][    ,     ]: choice set (1: t_1 = 35; 2: t_1 = 63)
#                   [[ ]][[x]][[ ]][    ,     ]: p-value combination
#                   [[ ]][[ ]][[x]][    ,     ]: a-value
#                   [[ ]][[ ]][[ ]][rows, cols]: matrix (rows: payment combination; cols: r-grid values)
# NOTE: p-value combinations
# p_0 = c(1.0, 1.0, 0.8, 0.5, 0.5, 0.4)
# p_1 = c(1.0, 0.8, 1.0, 0.5, 0.4, 0.5)
create_threshold_matrices = function(p_0, p_1, r_grid, omega, a_threshold_values, list_payments, t_0, t_1){
  
  list_threshold_mats = list()
  
  for (t_ in 1:length(t_1)){
    # Initialize matrix to store the thresholds
    mat_ = init_matrix_r_grid(list_payments, r_grid)
    
    # Initialize list to store the threshold matrices for all probability combinations
    list_threshold_mats[[t_]] = list()
    
    for (p_ in 1:length(p_0)){
      list_threshold_mats[[t_]][[p_]] = list()
      
      for (a_ in 1:length(a_threshold_values)){
        for (r_ in 1:length(r_grid)){
          # For given grid-point r, calculate the gamma-thresholds
          mat_[ , r_] = calc_gamma_threshold(p_0[p_], p_1[p_], r_grid[r_], 
                                                   omega, 
                                                   a_threshold_values[a_], 
                                                   list_payments[[t_]][[1]], 
                                                   list_payments[[t_]][[2]], 
                                                   t_0, t_1[t_])
        }  # end for r_
        
        # For given value of a, store the threshold matrix
        list_threshold_mats[[t_]][[p_]][[a_]] = mat_
        
      }  # end for a_
    }  # end for p_
  }  # end for t_
  
  return(list_threshold_mats)
}

# CREATE: Gamma threshold matrices
threshold_matrices =  create_threshold_matrices(p_0, p_1, r_grid, omega, a_threshold_values, list_payments, t_0, t_1)



# CALCULATION OF CHOICE PROBABILITIES -------------------------------------

# DEF FUNCTION: Create conditional choice probability matrices (for any given point in r_grid)
# Return: Matrices of conditional choice probabilities 
#         [[x]][[ ]][[ ]][    ,     ]: (t  , k) combination
#         [[ ]][[x]][[ ]][    ,     ]: (p_0, p_1) combination
#         [[ ]][[ ]][[x]][    ,     ]: specific a-value
#         [[ ]][[ ]][[x]][rows, cols]: matrix of conditional cdfs (row = payment, col = r-value)
create_cond_cdf_matrices = function(pars){
  
  # Initialize matrix to store the conditional choice probabilities
  cond_cdf_mats = threshold_matrices
  
  # Retrieve parameters
  mu_r        = pars[1]
  mu_gamma    = pars[2]
  sigma_r     = pars[3]
  sigma_gamma = pars[4]
  rho         = pars[5]
  
  # Calculate the conditional location and scale parameters
  mu_gamma_r    = mu_gamma_cond_r(r_grid, mu_r, mu_gamma, sigma_r, sigma_gamma, rho)
  sigma_gamma_r = sigma_gamma_cond_r(r_grid, sigma_gamma, rho)
  
  # Populate each cdf matrix
  for (t_ in 1:length(t_1)){
    for (p_ in 1:length(p_0)){
      for (a_ in 1:length(a_threshold_values)){
        for (r_ in 1:length(r_grid)){
          cond_cdf_mats[[t_]][[p_]][[a_]][ , r_] = pnorm(
            threshold_matrices[[t_]][[p_]][[a_]][ , r_], mu_gamma_r[r_], sigma_gamma_r[r_]
          )
          # cond_cdf_mats[[t_]][[p_]][[a_]][ , r_] = ptruncnorm(
          #   # threshold_matrices[[t_]][[p_]][[a_]][ , r_], c(-9, -9), c(9, 9), mu_gamma_r[r_], sigma_gamma_r[r_]
          #   threshold_matrices[[t_]][[p_]][[a_]][ , r_], c(-7, 0), c(5, 9), mu_gamma_r[r_], sigma_gamma_r[r_]
          # )
        }  # end for r_
      }  # end for a_
    }  # end for p_
  }  # end for t_
  
  return(cond_cdf_mats)
}


# DEF FUNCTION: Create weighted conditional choice probability matrices (for any given point in r_grid)
# Return: Matrices of weighted conditional choice probabilities 
#         [[x]][[ ]][[ ]][    ,     ]: (t  , k) combination
#         [[ ]][[x]][[ ]][    ,     ]: (p_0, p_1) combination
#         [[ ]][[ ]][[x]][    ,     ]: specific a-value
#         [[ ]][[ ]][[x]][rows, cols]: matrix of conditional cdfs (row = payment, col = r-value)
create_w_cond_cdf_matrices = function(pars){
  
  # Initialize matrix to store the conditional choice probabilities
  list_w_cond_cdf_mats = create_cond_cdf_matrices(pars)
  
  # Retrieve parameters
  mu_r        = pars[1]
  sigma_r     = pars[3]

  # Calculate the pdf-values for any grid point
  pdf_vals = dnorm(r_grid, mu_r, sigma_r)
  
  # Populate each weighted cdf matrix
  for (t_ in 1:length(t_1)){
    for (p_ in 1:length(p_0)){
      for (a_ in 1:length(a_threshold_values)){
        for (r_ in 1:length(r_grid)){
          list_w_cond_cdf_mats[[t_]][[p_]][[a_]][ , r_] =
            list_w_cond_cdf_mats[[t_]][[p_]][[a_]][ , r_] * pdf_vals[r_]
        }  # end for r_
      }  # end for a_
    }  # end for p_
  }  # end for t_
  
  return(list_w_cond_cdf_mats)
}


# DEF FUNCTION: Create the sum of weighted conditional choice probability matrices
# i.e., sum over the r-grid points
# Return: Matrices of choice probabilities
#         [[x]][[ ]][[ ]]: (t  , k) combination
#         [[ ]][[x]][[ ]]: (p_0, p_1) combination
#         [[ ]][[ ]][[x]]: specific a-value
create_sum_w_cdf_matrices = function(pars){
  
  # Create matrix of weighted choice probabilities (for any given r-grid point)
  mat_ = create_w_cond_cdf_matrices(pars)
  store_mat_ = list() 
  
  for (t_ in 1:length(t_1)){
    store_mat_[[t_]] = list()
    
    for (p_ in 1:length(p_0)){
      store_mat_[[t_]][[p_]] = matrix(NA, 
                                      nrow=length(list_payments[[1]][[1]]),
                                      ncol= length(a_threshold_values)
      )
      
      for (a_ in 1:length(a_threshold_values)){
        # Sum the weighted conditional choice probabilities across values of r_grid
        store_mat_[[t_]][[p_]][ , a_] = rowSums(mat_[[t_]][[p_]][[a_]])
        
        # Define sum-value for the special a-threshold a_S
        if (a_ == length(a_threshold_values)){
          store_mat_[[t_]][[p_]][ , a_] = calc_sum_pdf_vals(pars)
        }
        
      }  # end for a_
    }  # end for p_
  }  # end for t_
  
  return(store_mat_)
}


# DEF FUNCTION: Transform CDFs to probabilities
# Return: Probability matrices
# NOTE: This function needs to include the tremble probability
transform_cdf_to_pdfs_with_tremble = function(cdf_mats, pars){
  
  dims = dim(cdf_mats[[1]][[1]])
  # Retrieve tremble parameter
  tremble=pars[6]
  
  sum_pdf_vals = calc_sum_pdf_vals(pars)
  
  prob_mats = list()
  for (t_ in 1:length(t_1)){
    prob_mats[[t_]] = list()
    for (p_ in 1:length(p_0)){
      prob_mats[[t_]][[p_]] = matrix(
        NA, 
        nrow = dims[1], 
        ncol = dims[2]
      )
      
      for (as_index in 1:dims[2]){
        if(as_index == 1){
          prob_mats[[t_]][[p_]][, as_index] = 
            (1-tremble)*( (cdf_mats[[t_]][[p_]][ , as_index] - matrix(0, nrow=dims[1], ncol=1))/sum_pdf_vals ) + tremble/length(a_threshold_values)
        }  # end for as_index
        if (as_index > 1 & as_index <= length(a_threshold_values)){
          prob_mats[[t_]][[p_]][, as_index] = 
            (1-tremble)*( (cdf_mats[[t_]][[p_]][ , as_index] - cdf_mats[[t_]][[p_]][ , (as_index-1)])/sum_pdf_vals ) + tremble/length(a_threshold_values)
        } # end for as_index
        
      } # end for as_index
    } # end for p_
  } # end for t_
  
  return(prob_mats)
}



# DATA: LOAD EMPIRICAL DATA -----------------------------------------------

# Retrieve meta-information (individuals)
ids_individuals = readRDS(file = paste0(path_empirical_data, "id_list.rds"))
num_individuals = length(ids_individuals)


# DEF FUNCTION: Extract empirical data
# Return: List of empirical choice matrices for short and long horizon
get_empirical_data = function(path_data_short, path_data_long){
  
  readfile_short = readRDS(paste0(path_empirical_data, path_data_short, ".rds"))
  readfile_long  = readRDS(paste0(path_empirical_data, path_data_long, ".rds"))
  list_emp_data  = list(readfile_short, readfile_long)
  
  return(list_emp_data)
}

# Get empirical data
# NOTE: We do not remove the dominated options in the case of tremble
emp_data = get_empirical_data(
  "AS_empirical_matrices_1", 
  "AS_empirical_matrices_2"
)


# Initialize lists to store individual-level data
# NOTE: We do not remove the dominated options in the case of tremble
list_emp_data_individual = list()

for (id_ in seq_along(ids_individuals)){
  
  # Store individual-level data
  list_emp_data_individual[[id_]] = get_empirical_data(
    paste0("individual/AS_empirical_matrices_1_ind_", ids_individuals[id_]), 
    paste0("individual/AS_empirical_matrices_2_ind_", ids_individuals[id_]))
  
}



# ESTIMATION: DEFINE LOG-LIKELIHOOD FUNCTION ------------------------------

# DEF FUNCTION: Log-likelihood function with tremble
# Return: Value of the log-likelihood function
tremble_log_likelihood_function = function(pars, data){
  
  # Retrieve weighted CDF matrices
  sum_w_cdf_matrices = create_sum_w_cdf_matrices(pars)
  # Transform to choice probability matrices
  w_choice_probability_mats_ = transform_cdf_to_pdfs_with_tremble(sum_w_cdf_matrices, pars)

  # Initialize log-likelihood value
  likelihood = 0
  
  for (t_ in 1:length(t_1)){
    for (p_ in 1:length(p_0)){
      
      # Get log conditional choice probabilities
      w_cp_mats = w_choice_probability_mats_[[t_]][[p_]]
      log_w_cp_mats = log(w_cp_mats)
      log_cp = log_w_cp_mats
      
      # Set finite values whenever the log-likelihood function is not defined
      log_cp = apply(
        log_cp,
        c(1,2),
        (function(val) if(is.nan(val) | is.na(val)){val=-999999999} else{val})
      )
      # Set finite values whenever the log-likelihood function is infinity
      liks_ = apply(
        log_cp, 
        c(1,2), 
        (function(val) if(val==-Inf){val=-99999999999} else{val})
        )
      
      # Add up likelihood contributions
      liks = liks_ * data[[t_]][[p_]]
      sum_liks = sum(liks)
      
      # Set finite values whenever the sum of log-likelihood contributions is infinity
      if (sum_liks == - Inf){
        sum_liks = -999999999
      }
      
      # Add up likelihood contributions
      likelihood = likelihood + sum_liks
      
    }
  }
  
  # Set finite values whenever the sum of log-likelihood contributions is infinity (sanity check)
  if (likelihood == -Inf){
    likelihood = -999999999
  }
  
  print(paste0("Likelihood-value: ", likelihood))
  
  return(-likelihood)
}


    
# ESTIMATION: OPTIMIZATION ------------------------------------------------

# Initialize starting parameters
# NOTE: Starting parameter must include the tremble parmaeter now
starting_pars_ = c(0.1, 5, 1, 3, 0.7, 0.1)

# Set bounds
lower_bounds_  = c(-10, -10, 0.01, 0.01, -0.99, -1) 
upper_bounds_  = c( 10,  10,   10,   10,  0.99,  1)

# Optimization routine
mle = optim(
  par     = starting_pars_, 
  fn      = tremble_log_likelihood_function,
  lower   = lower_bounds_,
  upper   = upper_bounds_,
  data    = emp_data,
  method  = "L-BFGS-B", 
  hessian = TRUE
)



# ESTIMATION: CALCULATE STANDARD ERRORS -----------------------------------

# Extract gradients for each individual at the ML estimates and store them in a list 
list_mult_gradients_individual = list()
for (id_ in seq_along(ids_individuals)){
  
  # DEF FUNCTION: Wrapper for the log-likelihood function 
  # Return: Tremble log-likelihood function as function of the parameters only, given individual choices
  aux_tremble_log_likelihood_function = function(pars){
    
    aux_tremble_ll_function_ = tremble_log_likelihood_function(pars, list_emp_data_individual[[id_]])
    
    return(aux_tremble_ll_function_)
  }
  
  # Store the multiplied gradients in a list
  gradient                              = grad(aux_tremble_log_likelihood_function, mle$par)
  mult_gradients                        = gradient %*% t(gradient)
  list_mult_gradients_individual[[id_]] = mult_gradients
  
}

# Sum the gradients and get the inverse of the Hessian
sum_gradients_individual = Reduce("+", list_mult_gradients_individual)
inv_hessian              = solve(mle$hessian)

# Calculate robust standard errors, clustered at the individual level
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
  "mu_gamma"         = est_pars[2], "sigma_gamma"  = est_pars[4], 
  "rho"              = est_pars[5], "tremble"      = est_pars[6],
  "se_mu_r"          = est_se[1],  "se_sigma_r"    = est_se[3], 
  "se_mu_gamma"      = est_se[2], "se_sigma_gamma" = est_se[4], 
  "se_rho"           = est_se[5], 
  "likelihood_value" = mle$value, "omega"          = omega,
  "min_r_grid" = r_grid_info[1], "max_r_grid" = r_grid_info[2], "step_r_grid" = r_grid_info[3]
)

# Export parameter estimates
write.table(export_vector,
            file=paste0(path, "/Output/estimates/AS/", "CRRA_estimated_pars_omega_", omega, "_tremble", ".csv"),
            row.names=TRUE, col.names=FALSE, sep=",")
saveRDS(export_vector, file=paste0(path, "/Output/estimates/AS/", "CRRA_estimated_pars_omega_", omega, "_tremble", ".rds"))

