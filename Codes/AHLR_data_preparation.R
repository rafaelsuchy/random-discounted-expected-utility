# DESCRIPTION -------------------------------------------------------------

# This program prepares the AHLR data for estimation
# Reference data: Andersen et al. (2008)



# PATH AND PACKAGES -------------------------------------------------------

# Clean up environment
rm(list=ls())

# Load required packages
packages_required <- c("readxl", "MASS", "dplyr", "here")
install.packages( setdiff( packages_required, rownames( installed.packages() ) ) )
for(name in packages_required) { library(name, character.only = TRUE) }

# Set paths
path                = here::here()
path_empirical_data = paste0(path, "/Data/processed_data/AHLR/")

# Account for Indifferences?
# Account for (1) indifferences or discard (0) indifferences
account_indiffs = 1
# Set rate at which indifferences are accrued to the safe choice (lottery a)
rate_indiff = 0.5  



# EXTRACT CHOICE DATA -----------------------------------------------------

# Extract raw data
df_raw_data          = read_excel(paste0(path, "/Data/raw_data/AHLR_raw_data.xls"))
df_raw_data$DelayedP = round(df_raw_data$DelayedP)  # Transform DelayedDP to integer variable

# Read and save meta-data
id_list         = unique(df_raw_data$id)
num_individuals = length(id_list)
saveRDS(id_list, file = paste0(path_empirical_data, "id_list.rds"))


# COUNT RISK CHOICES ------------------------------------------------------

# DEF FUNCTION: Count choices in risk tasks
# Return: Empirical choice matrices for risk tasks
#         [[x]][    ,     ]: x="a"     - Number of choices for lottery a (safe)
#                            x="b"     - Number of choices for lottery b (risky)
#                            x="i"     - Number of indifferent choices (neither safe nor risky)
#                            x="total" - Total choices
#         [[ ]][rows, cols]: rows      - Probabilities for safe lottery 
#                            cols      - Payment combination
count_risk_choices = function(df_data){
  
  df_data_risk_task1 = df_data %>% 
    filter(Rtask1 == 1)
  df_data_risk_task2 = df_data %>%
    filter(Rtask2 == 1)
  df_data_risk_task3 = df_data %>% 
    filter(Rtask3 == 1)
  df_data_risk_task4 = df_data %>%
    filter(Rtask4 == 1)
  
  df_data_risk_tasks = list(
    df_data_risk_task1, df_data_risk_task2, 
    df_data_risk_task3, df_data_risk_task4
    )
  
  # Initialize empty matrices to store data
  empirical_matrix_risk_choice_a = matrix(
    data = NA,
    nrow = length(unique(na.omit(df_data$RowProbA))),
    ncol = length(df_data_risk_tasks)
  )
  empirical_matrix_risk_choice_b    = empirical_matrix_risk_choice_a
  empirical_matrix_risk_choice_i    = empirical_matrix_risk_choice_a
  empirical_matrix_risk_frequencies = empirical_matrix_risk_choice_a
  
  # Populate matrices of observed choices (a, b, indifferences)
  for (m_ in 1:length(df_data_risk_tasks)){
    
    risk_task = df_data_risk_tasks[[m_]]
    
    for (p_ in seq_along(unique(na.omit(df_data$RowProbA)))){
      p_risk_task = risk_task[round(risk_task$RowProbA, 2) == p_/10, ]
      choice_a = length(p_risk_task[p_risk_task$choice=="a", ]$choice)
      choice_b = length(p_risk_task[p_risk_task$choice=="b", ]$choice)
      choice_i = length(p_risk_task[p_risk_task$choice=="i", ]$choice)
      empirical_matrix_risk_choice_a[p_, m_] = choice_a + account_indiffs * rate_indiff * choice_i     # a: "safe"
      empirical_matrix_risk_choice_b[p_, m_] = choice_b + account_indiffs * (1-rate_indiff) * choice_i # b: "risky"
      empirical_matrix_risk_choice_i[p_, m_] = choice_i                                                # i: "indifference"
    }  # end for p_
  }  # end for m_
  
  # Count total choices
  total_risk_choices = empirical_matrix_risk_choice_a + empirical_matrix_risk_choice_b + 
    (1 - account_indiffs) * empirical_matrix_risk_choice_i

  # Collect empirical matrices
  empirical_matrices_risk = list(
    "a"=empirical_matrix_risk_choice_a, 
    "b"=empirical_matrix_risk_choice_b
    )
  
  return(empirical_matrices_risk)
}



# EXPORT RISK CHOICES -----------------------------------------------------

# Count risk choices
empirical_matrices_risk = count_risk_choices(df_raw_data)[c("a", "b")]

# Export data
write.table(
  as.data.frame(empirical_matrices_risk),
  file = paste0(path_empirical_data, "AHLR_empirical_matrices_risk", ".csv"),
  quote = F, sep = ",", row.names = F, col.names = TRUE)
saveRDS(empirical_matrices_risk, file = paste0(path_empirical_data, "AHLR_empirical_matrices_risk", ".rds"))



# COUNT TIME CHOICES ------------------------------------------------------

# DEF FUNCTION: Count number of choices in the time tasks
# Return: Empirical choice matrices for time tasks
#         [[x]][    ,     ]: x="a"     - Number of choices for earlier date
#                            x="b"     - Number of choices for later date
#                            x="i"     - Number of indifferent choices 
#                            x="total" - Total choices
#         [[ ]][rows, cols]: rows      - Return of later payment
#                            cols      - Time gap
count_time_choices = function(df_data){
  
  # Extract data for each time task
  df_data_time_horizon1 = df_data %>%
    filter(DRtask1 == 1)
  df_data_time_horizon4 = df_data %>%
    filter(DRtask2 == 1)
  df_data_time_horizon6 = df_data %>%
    filter(DRtask3 == 1)
  df_data_time_horizon12 = df_data %>%
    filter(DRtask4 == 1)
  df_data_time_horizon18 = df_data %>%
    filter(DRtask5 == 1)
  df_data_time_horizon24 = df_data %>%
    filter(DRtask6 == 1)
  
  df_data_time_menus = list(df_data_time_horizon1, df_data_time_horizon4, df_data_time_horizon6, 
                            df_data_time_horizon12, df_data_time_horizon18, df_data_time_horizon24)
  
  # Initialize empty matrices to store data
  empirical_matrix_time_choice_a = matrix(
    data = NA,
    nrow = length(unique(df_data_time_horizon1$DelayedP)),
    ncol = length(df_data_time_menus)
  )
  
  empirical_matrix_time_choice_b = empirical_matrix_time_choice_a
  empirical_matrix_time_choice_i = empirical_matrix_time_choice_a
  empirical_matrix_time_frequencies = empirical_matrix_time_choice_a
  
  # Populate matrices of observed choices
  for (m_ in seq_along(df_data_time_menus)){
    time_menu = df_data_time_menus[[m_]]
    cnt_ = 0
    
    for (p_ in sort(unique(time_menu$DelayedP))){
      cnt_ = cnt_ + 1
      p_time_menu = time_menu[time_menu$DelayedP == p_, ]
      choice_a = length(p_time_menu[p_time_menu$choice=="a", ]$choice)
      choice_b = length(p_time_menu[p_time_menu$choice=="b", ]$choice)
      choice_i = length(p_time_menu[p_time_menu$choice=="i", ]$choice)
      empirical_matrix_time_choice_a[cnt_, m_] = choice_a + account_indiffs * rate_indiff * choice_i     # a: "earlier"
      empirical_matrix_time_choice_b[cnt_, m_] = choice_b + account_indiffs * (1-rate_indiff) * choice_i # b: "later"
      empirical_matrix_time_choice_i[cnt_, m_] = choice_i                                                # i: "indifference"
    }  # end for p_
  }  # end for m_
  
  # Count total choices
  total_time_choices = empirical_matrix_time_choice_a + empirical_matrix_time_choice_b + (1 - account_indiffs) * empirical_matrix_time_choice_i
  
  # Collect empirical matrices for earlier (a) and later (b) choices
  empirical_matrices_time = list(
    "a" = empirical_matrix_time_choice_a, 
    "b" = empirical_matrix_time_choice_b
    )

  return(empirical_matrices_time)
}



# EXPORT TIME CHOICES -----------------------------------------------------

# Count time choices
empirical_matrices_time = count_time_choices(df_raw_data)[c("a", "b")]

# Export data
write.table(
  as.data.frame(empirical_matrices_time),
  file = paste0(path_empirical_data, "AHLR_empirical_matrices_time", ".csv"),
  quote=F, sep=",", row.names=F, col.names = TRUE
)
saveRDS(empirical_matrices_time, file = paste0(path_empirical_data, "AHLR_empirical_matrices_time", ".rds"))



# COUNT RISK AND TIME CHOICES BY INDIVIDUAL -------------------------------

# DEF FUNCTION: Filter for individuals by id
# Return: Dataframe that contains choices from a single subject
filter_individual = function(df_data, filter_id){
  
  filtered_df = df_data %>%
    filter(id == filter_id)
  
  return(filtered_df)
}

# Create output for each individual
for (id_ in 1:length(id_list)){
  
  # Extract raw data of each individual
  individual_data = filter_individual(df_raw_data, id_list[id_])
  
  # Count individual choices
  risk_choices_individual_data = count_risk_choices(individual_data)[c("a", "b")]
  time_choices_individual_data = count_time_choices(individual_data)[c("a", "b")]
  
  # Enlarge matrix (if smaller) with zero-choices
  risk_zero_choice_matrix = matrix(0, nrow = 4, ncol = ncol(risk_choices_individual_data$a))
  
  # Append zero-choices if the individual had less tasks than any other individual
  if (nrow(risk_choices_individual_data$a) < 9){
    
    risk_choices_individual_data$a = rbind(risk_choices_individual_data$a, risk_zero_choice_matrix) 
    risk_choices_individual_data$b = rbind(risk_choices_individual_data$b, risk_zero_choice_matrix)
  }
  
  # Save matrices
  # risk_filename_individual_ = paste0(path_empirical_data, "individual/", "AHLR_empirical_matrices_risk_ind_", id_list[id_])
  # time_filename_individual_ = paste0(path_empirical_data, "individual/", "AHLR_empirical_matrices_time_ind_", id_list[id_])

  # Export risk choices
  write.table(
    as.data.frame(risk_choices_individual_data),
    file = paste0(path_empirical_data, "individual/", "AHLR_empirical_matrices_risk_ind_", id_list[id_], ".csv"),
    quote = F, sep = ",", row.names = F, col.names = TRUE
  )
  saveRDS(risk_choices_individual_data, file = paste0(path_empirical_data, "individual/", "AHLR_empirical_matrices_risk_ind_", id_list[id_], ".rds"))
  
  # Export time choices
  write.table(
    as.data.frame(time_choices_individual_data),
    file = paste0(path_empirical_data, "individual/", "AHLR_empirical_matrices_time_ind_", id_list[id_], ".csv"),
    quote = F, sep = ",", row.names = F, col.names = TRUE
  )
  saveRDS(time_choices_individual_data, file = paste0(path_empirical_data, "individual/", "AHLR_empirical_matrices_time_ind_", id_list[id_], ".rds"))
}



# COUNT RISK AND TIME CHOICES BY CHARACTERISTICS --------------------------

# DEF FUNCTION: Filter dataframe for given characteristics
# Return: Dataframe that contains choices from subjects with given characteristics
filter_characteristics = function(df_data, filter_chars){

  filtered_df = df_data %>%
    filter(if ("young_middle" %in% filter_chars) young == 1 | middle == 1 else TRUE) %>%
    filter(if ("young_old" %in% filter_chars)    young == 1 | old == 1 else TRUE) %>%
    filter(if ("middle_old" %in% filter_chars)  middle == 1 | old == 1 else TRUE) %>%
    filter(if ("young"   %in% filter_chars)   young == 1 else TRUE) %>%
    filter(if ("semimiddle" %in% filter_chars) young == 0 & middle == 0 & old == 0 else TRUE) %>%
    filter(if ("middle"  %in% filter_chars)  middle == 1 else TRUE) %>%
    filter(if ("old"     %in% filter_chars)     old == 1 else TRUE) %>%
    filter(if ("female"  %in% filter_chars)  female == 1 else TRUE) %>%
    filter(if ("male"    %in% filter_chars) female == 0 else TRUE) %>%
    filter(if ("copen"   %in% filter_chars)   copen == 1 else TRUE) %>%
    filter(if ("city"    %in% filter_chars)    city == 1 else TRUE) %>%
    filter(if ("owner"   %in% filter_chars)   owner == 1 else TRUE) %>%
    filter(if ("retired" %in% filter_chars) retired == 1 else TRUE) %>%
    filter(if ("student" %in% filter_chars) student == 1 else TRUE) %>%
    filter(if ("skilled" %in% filter_chars) skilled == 1 else TRUE) %>%
    filter(if ("longedu" %in% filter_chars) longedu == 1 else TRUE) %>%
    filter(if ("kids"    %in% filter_chars)    kids == 1 else TRUE) %>%
    filter(if ("single"  %in% filter_chars)  single == 1 else TRUE) %>%
    filter(if ("IncLow"  %in% filter_chars)  IncLow == 1 else TRUE) %>%
    filter(if ("IncMed"  %in% filter_chars)  IncMed == 1 else TRUE) %>%
    filter(if ("IncHigh" %in% filter_chars) IncHigh == 1 else TRUE)
  
  return(filtered_df)
  
}

# CAUTION: It may happen that the procedure aborts because the characteristics combination is not available. 


# Define list of characteristics to be filtered for
# Complete list of characteristics:
#     female, young, semimiddle, middle, old, copen, city, owned, retired, student, skilled, longedu, kids, single, IncLow, IncMed, IncHigh

# List of characteristics (used in the paper)
list_chars = list(
  c("female"), c("male"),  # gender
  c("young"), c("semimiddle"), c("middle"), c("old")  # age groups
)

# Save characteristic lists
saveRDS(list_chars, file = paste0(path_empirical_data, "list_characteristics.rds"))

# Count choices by characteristics
overview_chars_id        = list()
overview_chars_id_num    = list()
number_individuals_chars = list()
empty_chars_combinations = list()
overview_chars_id_matrix = matrix(0, nrow=length(list_chars), ncol=4)
i_ = 0

for (chars_ in list_chars){
  i_ = i_ + 1
  
  # Filter data by characteristics
  chars_data = filter_characteristics(df_raw_data, chars_)
  
  # Count risk and time choices from filtered data 
  risk_choices_chars_data = count_risk_choices(chars_data)
  time_choices_chars_data = count_time_choices(chars_data)
  
  # Create variable for the characteristics - id mapping
  chars_name                          = paste(chars_, collapse = "_")
  chars_id_list                       = unique(chars_data$id)
  chars_id_number                     = length(chars_id_list)
  overview_chars_id[[chars_name]]     = paste(unique(chars_data$id), collapse = ", ")
  overview_chars_id_num[[chars_name]] = unique(chars_data$id)
  
  number_individuals_chars[[chars_name]] = chars_id_number
  
  if (nrow(chars_data) < 7){
    empty_chars_combinations[[chars_name]] = "Empty combination"
  }
  
  # Create filename for saving
  filename_risk_chars_ = paste0(path_empirical_data, "characteristics/", "AHLR_empirical_matrices_risk_chars_", chars_name)
  filename_time_chars_ = paste0(path_empirical_data, "characteristics/", "AHLR_empirical_matrices_time_chars_", chars_name)
  if (account_indiffs == 0){
    filename_risk_chars_ = paste0(filename_risk_chars_, "_wo_indiffs")
    filename_time_chars_ = paste0(filename_time_chars_, "_wo_indiffs")
  }
  
  # Export choices from risk tasks
  write.table(
    as.data.frame(risk_choices_chars_data),
    file = paste0(filename_risk_chars_, ".csv"),
    quote = F, sep = ",", row.names = F, col.names = TRUE
  )
  saveRDS(risk_choices_chars_data, file = paste0(filename_risk_chars_, ".rds"))
  
  # Export choices from time tasks
  write.table(
    as.data.frame(time_choices_chars_data),
    file = paste0(filename_time_chars_, ".csv"),
    quote = F, sep = ",", row.names = F, col.names = TRUE
  )
  saveRDS(time_choices_chars_data, file = paste0(filename_time_chars_, ".rds"))
  
  # Collect mapping "characteristics to id"
  overview_chars_id_matrix[i_, ] = c(names(overview_chars_id)[i_], chars_id_number, num_individuals, overview_chars_id[[i_]])
}


# Save list "number of individuals with each characteristics"
saveRDS(number_individuals_chars, file= paste0(path_empirical_data, "characteristics/", "number_individuals_characteristics.rds"))
write.table(
  as.data.frame(number_individuals_chars),
  file = paste0(path_empirical_data, "characteristics/", "number_individuals_characteristics.csv"),
  row.names=FALSE, col.names=TRUE, sep=",")

# Save overview matrix
write.table(
  as.data.frame(overview_chars_id_matrix),
  file = paste0(path_empirical_data, "characteristics/", "overview_characteristics_id_matrix.csv"),
  row.names = FALSE, col.names = FALSE, sep = ","
)

# Save overview mapping: characteristics to id
saveRDS(overview_chars_id_num, file= paste0(path_empirical_data, "characteristics/", "overview_characteristics_id.rds"))



# EXPLANATION CHARACTERISTICS ---------------------------------------------
# Explanations according to STATA codebook from raw data in Andersen et al. (2008)

# id: Subject id
# female: Being female (1: Female) [binary]
# young: Aged less than 30 (1: Less than 30) [binary]
# semimiddle: Aged between 30 and 40 (1: Between 30 and 40) [binary] - Note: This was added to capture all age ranges
# middle: Aged between 40 and 50 (1: Between 40 and 50) [binary]
# old: Aged over 50 (1: Over 50) [binary]
# copen: Lives in Copenhagen area (1: Copenhagen area) [binary]
# city: Lives in larger city of 20,000 or more (1: Lives in large cite) [binary]
# owner: Own home or apartment (1: Owns home or apartment) [binary]
# retired:  Retired (1: Retired) [binary]
# student: Student (1: Student) [binary]
# skilled: Some post-secondary education (1: Skilled) [binary]
# longedu: Substantial higher education	(1: Higher education) [binary]
# kids: Has children (1: Has children) [binary]
# single: Lives alone	(1: Living alone) [binary]
# IncLow: Lower level income	(1: Low income) [binary]
# IncMed: Medium level income	(1: Medium income) [binary]
# IncHigh: Higher level income (1: High income) [binary]