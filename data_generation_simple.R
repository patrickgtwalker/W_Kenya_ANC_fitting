generate_data_simple<-function(n_month,
                               n_sites,
                               survey_samples_site_month_min,
                               survey_samples_site_month_max,
                               ANC_samples_site_month_min,
                               ANC_samples_site_month_max,
                               pred_months
    
){

# Model parameters
survey_prev <- runif(n_month)  # Random uniform probabilities for each month
site_sd <- 2
log_odd_ratio_site <- rnorm(n_sites, 0, site_sd)  # Odds ratios for each site
ANC_log_odds_ratio<-rnorm(1,0,1)

# Initialize an empty dataframes for survey and ANC
survey_data <- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())
ANC_data<- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())

# Loop over sites and months
for (site in 1:n_sites) {
  for (month in 1:n_month) {
    survey_sample_size <- sample(survey_samples_site_month_min:survey_samples_site_month_max, 1)
    ANC_sample_size <- sample(ANC_samples_site_month_min:ANC_samples_site_month_max, 1)
    
    survey_prevalence_s_m<-odds_to_probability(probability_to_odds(survey_prev[month])*exp(log_odd_ratio_site[site]))
    ANC_prev_s_m<-odds_to_probability(probability_to_odds(survey_prevalence_s_m)*exp(ANC_log_odds_ratio))
    
    survey_positive <- rbinom(1, size = survey_sample_size, prob = survey_prevalence_s_m)
    ANC_samples<-rbinom(ANC_sample_size,1,ANC_prev_s_m)
    
    survey_data <- rbind(survey_data, data.frame(month=month, site=site,sample_size=survey_sample_size,positive=survey_positive))
    ANC_data <- rbind(ANC_data, data.frame(month=month, site=site, sample_size=1, positive=ANC_samples))
  }
}
ANC_data$ANC=1
survey_data$ANC=0

data_to_model<-rbind(ANC_data, survey_data)%>%
  filter(!(month %in% pred_months&ANC==0))

return(list(
  param_df=data.frame(preg_par=ANC_log_odds_ratio,
                  tau_site=site_sd),
  simulated_survey_data=survey_data,
  data_to_model=data_to_model
))
}
