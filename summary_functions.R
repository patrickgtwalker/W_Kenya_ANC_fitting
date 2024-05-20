

# Functions to convert between odds and probability
odds_to_probability <- function(odds) {
  # Calculate probability from given odds
  odds / (1 + odds)
}

probability_to_odds <- function(probability) {
  # Calculate odds from given probability
  probability / (1 - probability)
}


get_data_plots <- function(data) {
  cols<-viridis(2)
  # Extract prevalence by site
  prev_by_site <- data %>%
    group_by(site,ANC) %>%
    summarise(observed_mean = sum(positive) / sum(sample_size),
              observed_lower = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$lower,
              observed_upper = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$upper)
  
  prev_by_site_wide <- prev_by_site %>%
    pivot_wider(
      names_from = ANC,
      values_from = c(observed_mean, observed_lower, observed_upper),
      names_sep = "_",
      names_glue = "{.value}_{ifelse(ANC == 1, 'ANC', 'Survey')}"
    )

   # site plot setup
  site_plot <- ggplot(prev_by_site_wide, aes(x = observed_mean_Survey, y = observed_mean_ANC)) +
    geom_point() +
    geom_errorbar(aes(ymin = observed_lower_ANC, ymax = observed_upper_ANC), width = 0.01) +
    geom_errorbarh(aes(xmin = observed_lower_Survey, xmax = observed_upper_Survey), height = 0.01) +
    labs(x = "Survey Prevalence ", y = "ANC prevalence")+
    theme_minimal()+geom_abline()
  
  # Extract prevalence by month
   prev_by_month <- data %>%
     group_by(ANC,month) %>%
    summarise(observed_mean = sum(positive) / sum(sample_size),
              observed_lower = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$lower,
              observed_upper = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$upper)
   
   prev_by_month$Sample <- factor(prev_by_month$ANC, levels = c(0, 1), labels = c("Survey", "ANC"))
   
  
   temporal_plot <- ggplot(prev_by_month, aes(x = month,color=Sample)) +
     geom_point(aes(y = observed_mean)) +
     geom_errorbar(aes(ymin = observed_lower, ymax = observed_upper), width = 0.2)+
     labs(x = "Month", y = "Prevalence", title = "Observed Prevalence by Month") +
     theme_minimal() +
     xlim(0,max(data$month)) +
     scale_color_manual(
       name = "Sample", 
       values = c("Survey" = cols[1], "ANC" = cols[2]),  # Specify colors for "Survey" and "ANC"
       labels = c("Survey", "ANC")
     )
   
   
   # Composite plot
   overall_plot <- (temporal_plot | site_plot)  +
     plot_layout(widths = c(2, 1))
   
  return(list(
    site_plot = site_plot,
    temporal_plot = temporal_plot,
    overall_plot = overall_plot
  ))
  }


get_fit_summary <- function(fitted_model, survey_data, param_df=NULL, pred_months=NULL, pred_sites=NULL) {
  cols<-viridis(2)

  merged_test_validate<-fitted_model$data_w_probs%>%
    filter(ANC==0)%>%
    rename(train_sample_size=sample_size,
           train_positive=positive)%>%
    left_join(survey_data,by=c("month","site","ANC"))

  prev_by_month<-merged_test_validate%>%
    select(month,site,sample_size,positive,starts_with("V"))%>%
    pivot_longer(cols = starts_with("V"), names_to = "prediction", values_to = "value")%>%
    group_by(month, prediction) %>%
    summarise(weighted_sum = sum(value * sample_size), 
              total_positive = sum(positive), 
              total_sample_size = sum(sample_size), 
              weighted_prediction = weighted_sum / total_sample_size,
              observed_proportion=total_positive/total_sample_size)%>%
    summarise( model_lower = quantile(weighted_prediction, probs = 0.025),
               model_median = quantile(weighted_prediction, probs = 0.5),
               model_upper = quantile(weighted_prediction, probs = 0.975),
               total_positive = first(total_positive),
               total_sampled = first(total_sample_size))%>%
    mutate(observed_mean = total_positive / total_sampled,
           observed_lower = binom.confint(total_positive,total_sampled, methods = "wilson")$lower,
           observed_upper = binom.confint(total_positive, total_sampled, methods = "wilson")$upper
    )%>%
    mutate(Status = ifelse(month %in% pred_months, "Predicted", "Fitted")) %>%
    arrange(month) %>%
    mutate(change_in_month = c(1, diff(month) > 1),
           change_in_status = c(1, diff(as.numeric(as.factor(Status))) != 0),
           Segment = cumsum(change_in_month | change_in_status)) %>%
    select(-change_in_month, -change_in_status) 
   
  ## plot comparison
  temporal_plot <- ggplot(prev_by_month, aes(x = month)) +
    geom_point(aes(y = observed_mean), color = "blue") +
    geom_errorbar(aes(ymin = observed_lower, ymax = observed_upper), width = 0.2, color = "blue") +
    geom_ribbon(aes(ymin = model_lower, ymax = model_upper, fill = Status, group = Segment), alpha = 0.5) +
    geom_line(aes(y = model_median, color = Status, group = Segment)) +
    scale_fill_manual(values = c("Predicted" = cols[1], "Fitted" = cols[2])) +
    scale_color_manual(values = c("Predicted" = cols[1], "Fitted" = cols[2])) +
    labs(x = "Month", y = "Prevalence", title = "Observed vs Modelled Prevalence by Month") +
    theme_minimal() +
    guides(fill = guide_legend(title = "Legend"), color = guide_legend(title = "Legend"))
  
  ### get observed (fitted and predicted) prevalence by site
  prev_by_site<-merged_test_validate%>%
    select(month,site,sample_size,positive,starts_with("V"))%>%
    pivot_longer(cols = starts_with("V"), names_to = "prediction", values_to = "value")%>%
    group_by(site, prediction) %>%
    summarise(weighted_sum = sum(value * sample_size), 
              total_positive = sum(positive), 
              total_sample_size = sum(sample_size), 
              weighted_prediction = weighted_sum / total_sample_size,
              observed_proportion=total_positive/total_sample_size)%>%
    summarise( model_lower = quantile(weighted_prediction, probs = 0.025),
               model_median = quantile(weighted_prediction, probs = 0.5),
               model_upper = quantile(weighted_prediction, probs = 0.975),
               total_positive = first(total_positive),
               total_sampled = first(total_sample_size))%>%
    mutate(observed_mean = total_positive / total_sampled,
           observed_lower = binom.confint(total_positive,total_sampled, methods = "wilson")$lower,
           observed_upper = binom.confint(total_positive, total_sampled, methods = "wilson")$upper,
           Status = ifelse(site %in% pred_sites, "Predicted", "Fitted")
           )
  
  
  ##make plot
  site_plot <- ggplot(prev_by_site, aes(x = observed_mean, y = model_median,color=Status)) +
    geom_point() +
    geom_errorbar(aes(ymin = model_lower, ymax = model_upper), width = 0.01) +
    geom_errorbarh(aes(xmin = observed_lower, xmax = observed_upper), height = 0.01) +
    labs(x = "Observed Prevalence (Mean)", y = "Modeled Prevalence (Median)")+
    theme_minimal()+geom_abline()
  
  
  # Posterior distribution analysis
  sample_matrix <- as.matrix(fitted_model$draws)
  post_samples <- as.data.frame(sample_matrix[, !grepl("^theta", colnames(sample_matrix))])
  posterior_long <- pivot_longer(post_samples, cols = everything(), names_to = "parameter", values_to = "value")
  color_mapping <- viridis::viridis(length(unique(posterior_long$parameter)), option = "C")
  names(color_mapping) <- unique(posterior_long$parameter)
  posterior_long$color <- color_mapping[as.character(posterior_long$parameter)]
  
  if(is.null(param_df)){
    posterior_plot <- ggplot(posterior_long, aes(x = value, y = parameter, fill = color)) +
      geom_density_ridges(rel_min_height = 0.01, alpha = 0.25) +
      scale_fill_identity() +
      scale_color_identity() +
      theme_ridges() +
      labs(title = "Posterior Distributions with Simulated Values", x = "Parameter Values", y = "Parameter") +
      theme(legend.position = "none")
  }else{
  simulated_long <- pivot_longer(param_df, cols = everything(), names_to = "parameter", values_to = "value")
  simulated_long$color <- color_mapping[as.character(simulated_long$parameter)]
 
  posterior_plot <- ggplot(posterior_long, aes(x = value, y = parameter, fill = color)) +
    geom_density_ridges(rel_min_height = 0.01, alpha = 0.25) +
    scale_fill_identity() +
    geom_vline(data = simulated_long, aes(xintercept = value, color = color), linetype = "dashed", linewidth = 1) +
    scale_color_identity() +
    theme_ridges() +
    labs(title = "Posterior Distributions with Simulated Values", x = "Parameter Values", y = "Parameter") +
    theme(legend.position = "none")
  }
  # Composite plot
  overall_plot <- (posterior_plot | site_plot) / temporal_plot +
    plot_layout(heights = c(1, 2))
  
  if(is.null(pred_months)&is.null(pred_sites)){
  return(list(
    posterior_plot = posterior_plot,
    site_plot = site_plot,
    temporal_plot = temporal_plot,
    overall_plot = overall_plot
  ))
  }
  else{
    data_to_predict<-merged_test_validate%>%filter(month%in%pred_months|site%in%pred_sites)
    prediction_metrics<-get_metrics(data_to_predict)
    return(list(
      posterior_plot = posterior_plot,
      site_plot = site_plot,
      temporal_plot = temporal_plot,
      overall_plot = overall_plot,
      prediction_metrics=prediction_metrics
    ))
    }
}

get_metrics<-function(data_to_predict){

  ### calculate site-month specific CRPS
  site_month_summary <- data_to_predict %>%
    rowwise() %>%
    mutate(
      CRPS = mean(crps_binom(
        y = positive, 
        size = sample_size, 
        prob = c_across(starts_with("V"))
      )),
      median=quantile(c_across(starts_with("V")),0.5),
      point_prev=positive/sample_size,
      RMSE=((median-point_prev)^2)^0.5
    ) %>%
    ungroup() %>%
    select(site,month,CRPS,median,point_prev,RMSE)
  
  month_summary <- data_to_predict %>%
    mutate(across(starts_with("V"), ~ . * sample_size)) %>%
    group_by(month) %>%
    summarise(
      across(starts_with("V"), sum, .names = "sum_{.col}"),
      sample_size = sum(sample_size),
      positive = sum(positive)
    ) %>%
    mutate(across(starts_with("sum_V"), ~ ./sample_size, .names = "weighted_{.col}"))%>%
    rowwise() %>%
    mutate(
      site="all",
      CRPS = mean(crps_binom(
        y = positive, 
        size = sample_size, 
        prob = c_across(starts_with("weighted_sum_V"))
      )),
      median=quantile(c_across(starts_with("weighted_sum_V")),0.5),
      point_prev=positive/sample_size,
      RMSE=((median-point_prev)^2)^0.5
    ) %>%
    ungroup() %>%
    select(site,month,CRPS,median,point_prev,RMSE)

  site_summary<-data_to_predict %>%
    mutate(across(starts_with("V"), ~ . * sample_size)) %>%
    group_by(site) %>%
    summarise(
      across(starts_with("V"), sum, .names = "sum_{.col}"),
      sample_size = sum(sample_size),
      positive = sum(positive)
    ) %>%
    mutate(across(starts_with("sum_V"), ~ ./sample_size, .names = "weighted_{.col}"))%>%
    rowwise() %>%
    mutate(
      month="all",
      CRPS = mean(crps_binom(
        y = positive, 
        size = sample_size, 
        prob = c_across(starts_with("weighted_sum_V"))
      )),
      median=quantile(c_across(starts_with("weighted_sum_V")),0.5),
      point_prev=positive/sample_size,
      RMSE=((median-point_prev)^2)^0.5
    ) %>%
    ungroup() %>%
    select(site,month,CRPS,median,point_prev,RMSE)

  summary_pred<-rbind(site_month_summary,month_summary,site_summary)
  
  return(
    summary_pred
  )
}