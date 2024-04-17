'%notin%' <- function(x, table) {
  !(x %in% table)
}

# Functions to convert between odds and probability
odds_to_probability <- function(odds) {
  # Calculate probability from given odds
  odds / (1 + odds)
}

probability_to_odds <- function(probability) {
  # Calculate odds from given probability
  probability / (1 - probability)
}

# Function to calculate weighted averages of predictions across months
get_weight_month_pred <- function(fit_array, pop_survey_data) {
  # Determine the number of months and predictions
  n_month <- dim(fit_array)[3]
  n_site <- dim(fit_array)[2]
  n_preds <- dim(fit_array)[1]
  # Initialize an array for storing weighted averages
  weighted_averages <- array(dim = c(n_month, n_preds))
  
  # Create a grid of all possible site-month combinations
  full_grid <- expand.grid(site = 1:n_site, month = 1:n_month)
  
  # Merge this grid with population survey data to include sample sizes
  # Replace missing values with 0, indicating no samples
  sample_sizes_df <- full_grid %>%
    left_join(pop_survey_data, by = c("site", "month")) %>%
    replace_na(list(sample_size = 0)) %>%
    arrange(site, month)
  
  # Convert the merged data into a matrix format for easier processing
  sample_sizes_matrix <- matrix(sample_sizes_df$sample_size, nrow = n_site, byrow = TRUE)
  
  # Loop through each prediction and month to calculate weighted averages
  for (pred in 1:n_preds) {
    for (month in 1:n_month) {
      predictions <- fit_array[pred, , month]  # Extract predictions for all sites
      weights <- sample_sizes_matrix[, month]  # Extract corresponding weights
      # Calculate the weighted average and store it
      weighted_averages[month, pred] <- sum(predictions * weights) / sum(weights)
    }
  }
  return(weighted_averages)
}

# Function to calculate weighted averages of predictions across sites
get_weight_site_pred <- function(fit_array, pop_survey_data) {
  # The dimensions for n_month should be corrected to n_site for consistency in variable naming
  n_month <- dim(fit_array)[3]
  n_site <- dim(fit_array)[2]
  n_preds <- dim(fit_array)[1]
  # Initialize an array for storing weighted averages
  weighted_averages <- array(dim = c(n_site, n_preds))
  
  # Create a grid of all possible site-month combinations
  full_grid <- expand.grid(site = 1:n_site, month = 1:n_month)
  
  # Merge this grid with population survey data to include sample sizes
  sample_sizes_df <- full_grid %>%
    left_join(pop_survey_data, by = c("site", "month")) %>%
    replace_na(list(sample_size = 0)) %>%
    arrange(site, month)
  
  # Convert the dataframe to a matrix of sample sizes
  sample_sizes_matrix <- matrix(sample_sizes_df$sample_size, nrow = n_site, byrow = TRUE)
  
  # Loop through each prediction and site to calculate weighted averages
  for (pred in 1:n_preds) {
    for (site in 1:n_site) {
      predictions <- fit_array[pred, site, ]  # Extract predictions for all months
      weights <- sample_sizes_matrix[site, ]  # Extract corresponding weights
      # Calculate the weighted average and store it
      weighted_averages[site, pred] <- sum(predictions * weights) / sum(weights)
    }
  }
  return(weighted_averages)
}



