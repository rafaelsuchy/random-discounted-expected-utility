# DESCRIPTION -------------------------------------------------------------

# This program prepares the AS data for estimation
# Reference data: Andreoni and Sprenger (2012)



# PATH AND PACKAGES -------------------------------------------------------

# Clean up environment
rm(list=ls())

# Load required packages
packages_required <- c("MASS", "pracma", "msm", "here", "readxl", "dplyr")
install.packages( setdiff( packages_required, rownames( installed.packages() ) ) )
for(name in packages_required) { library(name, character.only = TRUE) }

# Set paths
path                = here::here()
path_empirical_data = paste0(path, "/Data/processed_data/AS/")



# EXPERIMENTAL PARAMETERS -------------------------------------------------

# Probability sequences
p_0 = c(1.0, 1.0, 0.8, 0.5, 0.5, 0.4)
p_1 = c(1.0, 0.8, 1.0, 0.5, 0.4, 0.5)

# Timings
t_0 = c(7)           # after 7 days
t_1 = c(28+7, 56+7)  # after 35 / 63 days 

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

# Define a-thresholds
a_thresholds = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1)



# EXTRACT CHOICE DATA -----------------------------------------------------

# Extract raw data, create column that track budget allocation a (a = allocation to t_1)
df_data_raw     = read_excel(paste0(path, "/Data/raw_data/AS_raw_data.xls"))
df_data_raw$v_a = (100 - df_data_raw$v) / 100  

# Select relevant variables only
df_data = df_data_raw[c(
  "labnumber", "budgetnum", "k", "perc1", "perc2", "rate1", "cons1", "cons2", 
  "v", "v_a", "risk100100", "risk10080", "risk5050", "risk5040", "risk4050"
  )] 

# Subset dataframe into short and long horizons
list_emp_data = list(
  df_data_28 = df_data[df_data$k == "k = 28 days", ],
  df_data_56 = df_data[df_data$k == "k = 56 days", ]
  )

# Read and save meta-data
id_list = unique(df_data$labnumber)
saveRDS(id_list, file=paste0(path, "/Data/processed_data/AS/", "id_list.rds"))



# COUNT CHOICES -----------------------------------------------------------

# CREATE: List of empirical choice matrices
emp_data_a_thresholds = list()
for (t_ in 1:length(t_1)){
  emp_data_a_thresholds[[t_]] = list()
  df_emp_data_ = list_emp_data[[t_]]
  
  for (p_ in 1:length(p_0)){
    selector_p_ = (df_emp_data_$perc1/100) == p_0[p_] & (df_emp_data_$perc2/100) == p_1[p_] 
    df_emp_data_cond_p_ = df_emp_data_[selector_p_, ]
    
    mat_ = matrix(NA, nrow=length(list_payments[[1]][[1]]), ncol=length(a_thresholds))
    
    for (a_ in 1:length(a_thresholds)){
      
      # Assign "a=0" choices
      if (a_ == 1){ 
        mat_[, a_] = c(
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 20 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 19 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 18 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 17 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 16 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 15 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 14 & df_emp_data_cond_p_$v_a < a_thresholds[a_])
        )
      }  
      
      # Asssign "interior choices"
      if (a_ > 1 & a_ < length(a_thresholds)){
        mat_[, a_] = c(
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 20 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 19 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 18 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 17 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 16 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 15 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 14 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_])
        )
      }  
      
      if (a_ == length(a_thresholds)){  # Assign a=S choices
        mat_[, a_] = c(
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 20 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 19 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 18 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 17 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 16 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 15 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
          sum(round(df_emp_data_cond_p_$rate1*100, 0) == 14 & df_emp_data_cond_p_$v_a >= a_thresholds[a_])
        )
      } 
      
    } # end for a_
    
    emp_data_a_thresholds[[t_]][[p_]] = mat_
    
  }  # end for p_
}  # end for t_

# Structure of "emp_data_a_thresholds"
# [[x]][[ ]]: index for (t_0, t_1) combination
# [[ ]][[x]]: index for (p_0, p_1) combination

# Available (p_0, p_1) combinations
# 1 = (1.0, 1.0); 
# 2 = (1.0, 0.8); 
# 3 = (0.8, 1.0); 
# 4 = (0.5, 0.5); 
# 5 = (0.5, 0.4); 
# 6 = (0.4, 0.5)


# Save empirical choice matrix
for (file_ in 1:length(emp_data_a_thresholds)){
  write.table(
    as.data.frame(emp_data_a_thresholds[[file_]]),
    file=paste0(path, "/Data/processed_data/AS/", "AS_empirical_matrices_", file_, ".csv"),
    quote=F, sep=",",row.names=FALSE, col.names = TRUE)
  saveRDS(
    emp_data_a_thresholds[[file_]], 
    file=paste0(path, "/Data/processed_data/AS/", "AS_empirical_matrices_", file_, ".rds"))
}



# COUNT CHOICES BY INDIVIDUAL ---------------------------------------------

# DEF FUNCTION: Filter for individuals by id
# Return: Dataframe that contains choices from a single subject
filter_individual = function(df_data, filter_id){
  
  filtered_df = df_data %>%
    filter(labnumber == filter_id)
  
  return (filtered_df)
  
}


# Count data by individual
for (id in id_list){
  
  list_emp_data_individual = list(
    filter_individual(list_emp_data$df_data_28, id), 
    filter_individual(list_emp_data$df_data_56, id)
  )
  
  list_emp_data_a_thresholds = list()
  for (t_ in 1:length(t_1)){
    list_emp_data_a_thresholds[[t_]] = list()
    df_emp_data_ = list_emp_data_individual[[t_]]

    for (p_ in 1:length(p_0)){
      
      # Subset data with respect to probability combination
      selector_p_ = (df_emp_data_$perc1/100) == p_0[p_] & (df_emp_data_$perc2/100) == p_1[p_] 
      df_emp_data_cond_p_ = df_emp_data_[selector_p_ , ]
      
      # Initialize matrix for choices by thresholds
      mat_ = matrix(NA, nrow=length(list_payments[[1]][[1]]), ncol=length(a_thresholds))
      
      for (a_ in 1:length(a_thresholds)){
        
        # Assign "a=0" choices
        if (a_ == 1){ 
          mat_[, a_] = c(
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 20 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 19 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 18 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 17 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 16 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 15 & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 14 & df_emp_data_cond_p_$v_a < a_thresholds[a_])
          )
        }  # end if
        
        # Assign "interior choices"
        if (a_ > 1 & a_ < length(a_thresholds)){
          mat_[, a_] = c(
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 20 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 19 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 18 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 17 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 16 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 15 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 14 & df_emp_data_cond_p_$v_a >= a_thresholds[(a_-1)] & df_emp_data_cond_p_$v_a < a_thresholds[a_])
          )
        }  # end if
        
        # Assign "a=S" choices
        if (a_ == length(a_thresholds)){ 
          mat_[, a_] = c(
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 20 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 19 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 18 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 17 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 16 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 15 & df_emp_data_cond_p_$v_a >= a_thresholds[a_]),
            sum(round(df_emp_data_cond_p_$rate1*100, 0) == 14 & df_emp_data_cond_p_$v_a >= a_thresholds[a_])
          )
        } # end if
        
      } # end for a_
      
      list_emp_data_a_thresholds[[t_]][[p_]] = mat_
      
    }  # end for p_
  }  # end for t_
  
  for (file_ in 1: length(list_emp_data_a_thresholds)){
    write.table(
      as.data.frame(list_emp_data_a_thresholds[[file_]]),
      file=paste0(path, "/Data/processed_data/AS/individual/", "AS_empirical_matrices_", file_, "_ind_", id, ".csv"),
      quote=F, sep=",",row.names=FALSE, col.names = TRUE)
    saveRDS(
      list_emp_data_a_thresholds[[file_]], 
      file=paste0(path, "/Data/processed_data/AS/individual/", "AS_empirical_matrices_", file_, "_ind_", id, ".rds")
    )
  }  # end for file_
}  # end for id