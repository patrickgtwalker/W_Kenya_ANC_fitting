generate_data_m1<-function(n_month=50,
                                   n_sites=20,
                               base_par=0,
                               month_sd=0.6,
                               site_sd =0.8,
                                   ANC_log_odds_ratio=-0.5,
                                   survey_samples_site_month_min=20,
                                   survey_samples_site_month_max=100,
                                   ANC_samples_site_month_min=10,
                                   ANC_samples_site_month_max=50
                                   
){
  
  # generate prevalence and site-specific params
  log_odd_ratio_month <- rnorm(n_month, 0, month_sd)  # Odds ratios for each site
  log_odd_ratio_site <- rnorm(n_sites, 0, site_sd)  # Odds ratios for each site

  # Initialize an empty dataframes for survey and ANC
  survey_data <- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())
  ANC_data<- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())
  
  # Loop over sites and months
  for (site in 1:n_sites) {
    for (month in 1:n_month) {
      survey_sample_size <- sample(survey_samples_site_month_min:survey_samples_site_month_max, 1)
      ANC_sample_size <- sample(ANC_samples_site_month_min:ANC_samples_site_month_max, 1)
      odds_survey<-exp(base_par+log_odd_ratio_month[month]+log_odd_ratio_site[site])
      survey_prevalence_s_m<-odds_to_probability(odds_survey)
      odds_ANC<-exp(base_par+log_odd_ratio_month[month]+log_odd_ratio_site[site]+ANC_log_odds_ratio)
      ANC_prev_s_m<-odds_to_probability(odds_ANC)
      survey_positive <- rbinom(1, size = survey_sample_size, prob = survey_prevalence_s_m)
      ANC_samples<-rbinom(ANC_sample_size,1,ANC_prev_s_m)
      
      survey_data <- rbind(survey_data, data.frame(month=month, site=site,sample_size=survey_sample_size,positive=survey_positive))
      ANC_data <- rbind(ANC_data, data.frame(month=month, site=site, sample_size=1, positive=ANC_samples))
    }
  }
  ANC_data$ANC=1
  survey_data$ANC=0
  
  data_to_model<-rbind(ANC_data, survey_data)
  
  return(list(
    param_df=data.frame(
      base_par=base_par,
      preg_par=ANC_log_odds_ratio,
      tau_month=month_sd,
      tau_site=site_sd),
    simulated_survey_data=survey_data,
    data_to_model=data_to_model
  ))
}

generate_data_m1_runif_site<-function(n_month=50,
                                       n_sites=20,
                                       month_sd =1,
                                       ANC_log_odds_ratio=-0.5,
                                       survey_samples_site_month_min=20,
                                       survey_samples_site_month_max=100,
                                       ANC_samples_site_month_min=10,
                                       ANC_samples_site_month_max=50,
                                       pred_sites=c(15:20)
                                       
){
  
  # generate prevalence and site-specific params
  survey_prev <- runif(n_sites)  # Random uniform probabilities for each site
  log_odd_ratio_month <- rnorm(n_month, 0, month_sd)  # Odds ratios for each site
  
  
  # Initialize an empty dataframes for survey and ANC
  survey_data <- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())
  ANC_data<- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())
  
  # Loop over sites and months
  for (site in 1:n_sites) {
    for (month in 1:n_month) {
      survey_sample_size <- sample(survey_samples_site_month_min:survey_samples_site_month_max, 1)
      ANC_sample_size <- sample(ANC_samples_site_month_min:ANC_samples_site_month_max, 1)
      
      survey_prevalence_s_m<-odds_to_probability(probability_to_odds(survey_prev[site])*exp(log_odd_ratio_month[month]))
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
    mutate(
      sample_size = ifelse(ANC == 0 & site %in% pred_sites, 0, sample_size),
      positive = ifelse(ANC == 0 & site %in% pred_sites, 0, positive)
    )
  
  
  
  return(list(
    param_df=data.frame(preg_par=ANC_log_odds_ratio,
                        tau_month=month_sd),
    pred_sites=pred_sites,
    simulated_survey_data=survey_data,
    data_to_model=data_to_model
  ))
}

generate_data_m1_runif_month<-function(n_month=50,
                               n_sites=20,
                               site_sd =1,
                               ANC_log_odds_ratio=-0.5,
                               survey_samples_site_month_min=20,
                               survey_samples_site_month_max=100,
                               ANC_samples_site_month_min=10,
                               ANC_samples_site_month_max=50,
                               pred_months=c(4:28,32:50)
    
){

# generate prevalence and site-specific params
survey_prev <- runif(n_month)  # Random uniform probabilities for each month
log_odd_ratio_site <- rnorm(n_sites, 0, site_sd)  # Odds ratios for each site


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
  mutate(
    sample_size = ifelse(ANC == 0 & month %in% pred_months, 0, sample_size),
    positive = ifelse(ANC == 0 & month %in% pred_months, 0, positive)
  )

return(list(
  param_df=data.frame(preg_par=ANC_log_odds_ratio,
                  tau_site=site_sd),
  pred_months=pred_months,
  simulated_survey_data=survey_data,
  data_to_model=data_to_model
))
}


generate_data_m5<-function(n_month=50,
                               n_sites=20,
                               base_par=0,
                               month_sd=0.8,
                               site_sd =1,
                               ANC_log_odds_intercept_cat1=0.1,
                               ANC_log_odds_intercept_cat2=-0.2,
                               ANC_log_odds_intercept_cat3=-0.5,
                               ANC_log_odds_gradient_cat1=-0.15,
                               ANC_log_odds_gradient_cat2=-0.3,
                               ANC_log_odds_gradient_cat3=-0.6,
                               prop_cat1=0.2,
                               prop_cat2=0.35,
                               survey_samples_site_month_min=20,
                               survey_samples_site_month_max=100,
                               ANC_samples_site_month_min=10,
                               ANC_samples_site_month_max=50
                               
){
  
  # generate prevalence and site-specific params
  #survey_prev <- runif(n_month)  # Random uniform probabilities for each month
  log_odd_ratio_month <- rnorm(n_month, 0, month_sd)  # Odds ratios for each site
  log_odd_ratio_site <- rnorm(n_sites, 0, site_sd)  # Odds ratios for each site
  
  
  # Initialize an empty dataframes for survey and ANC
  survey_data <- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())
  ANC_data<- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())
  
  # Loop over sites and months
  for (site in 1:n_sites) {
    for (month in 1:n_month) {
      survey_sample_size <- sample(survey_samples_site_month_min:survey_samples_site_month_max, 1)
      ANC_sample_size <- sample(ANC_samples_site_month_min:ANC_samples_site_month_max, 1)
      ANC_sample_by_grav_cat<-rmultinom(ANC_sample_size,1,prob=c(prop_cat1,prop_cat2,1-prop_cat1-prop_cat2))
      gravcat1<-ANC_sample_by_grav_cat[1,]
      gravcat2<-ANC_sample_by_grav_cat[2,]
      gravcat3<-ANC_sample_by_grav_cat[3,]
      #p<-ilogit(base_par+primi_trans_par*theta_site[as.numeric(adjust_data$site)]*adjust_data$primi+primi_par*adjust_data$primi+multi_trans_par*(theta_site[as.numeric(adjust_data$site)])*adjust_data$multi+multi_par*adjust_data$multi+grand_trans_par*(theta_site[as.numeric(adjust_data$site)])*adjust_data$grand+grand_par*adjust_data$grand+theta_month[as.numeric(adjust_data$month_code)]+theta_site[as.numeric(adjust_data$site)])
      survey_prevalence_s_m<-odds_to_probability(exp(base_par+log_odd_ratio_month[month]+log_odd_ratio_site[site]))
      odds_ANC<-exp(
        base_par+log_odd_ratio_month[month]+log_odd_ratio_site[site]+(ANC_log_odds_intercept_cat1+log_odd_ratio_site[site]*ANC_log_odds_gradient_cat1)*gravcat1+
        (ANC_log_odds_intercept_cat2+log_odd_ratio_site[site]*ANC_log_odds_gradient_cat2)*gravcat2+
        (ANC_log_odds_intercept_cat3+log_odd_ratio_site[site]*ANC_log_odds_gradient_cat3)*gravcat3
      )
      ANC_prev_s_m<-odds_to_probability(odds_ANC)
      survey_positive <- rbinom(1, size = survey_sample_size, prob = survey_prevalence_s_m)
      ANC_samples<-rbinom(n=ANC_sample_size,size=1,prob=ANC_prev_s_m)
      
      survey_data <- rbind(survey_data, data.frame(month=month, site=site,sample_size=survey_sample_size,positive=survey_positive,gravcat1=0,gravcat2=0,gravcat3=0))
      ANC_data <- rbind(ANC_data, data.frame(month=month, site=site, sample_size=1, positive=ANC_samples,gravcat1=gravcat1,gravcat2=gravcat2,gravcat3=gravcat3))
    }
  }
  ANC_data$ANC=1
  survey_data$ANC=0
  
  data_to_model<-rbind(ANC_data, survey_data)
  
  return(list(
    param_df=data.frame(
      base_par=base_par,
      tau_site=site_sd,
      tau_month=month_sd,
      log_odd_ratio_intercept_GC1=ANC_log_odds_intercept_cat1,
                        log_odd_ratio_intercept_GC2=ANC_log_odds_intercept_cat2,
                        log_odd_ratio_intercept_GC3=ANC_log_odds_intercept_cat3,
                        log_odd_ratio_gradient_GC1=ANC_log_odds_gradient_cat1,
                        log_odd_ratio_gradient_GC2=ANC_log_odds_gradient_cat2,
                        log_odd_ratio_gradient_GC3=ANC_log_odds_gradient_cat3
                        ),
    simulated_survey_data=survey_data,
    data_to_model=data_to_model
  ))
}

