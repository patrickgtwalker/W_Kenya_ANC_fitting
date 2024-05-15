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


get_fit_plots <- function(fitted_model, survey_data, param_df, pred_months=F, pred_sites=F) {
 
  cols<-viridis(2)
  
  # Extract predictions and observed statistics
  prev_by_site <- survey_data %>%
    group_by(site) %>%
    summarise(observed_mean = sum(positive) / sum(sample_size),
              observed_lower = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$lower,
              observed_upper = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$upper)%>%
    mutate(Status = ifelse(site %in% pred_sites, "Predicted", "Fitted"))
  
  fit_by_site<-fitted_model$data_w_probs%>%
    filter(ANC==0)%>%
    tidyr::pivot_longer(cols = starts_with("V"), names_to = "variable", values_to = "prediction") %>%
    mutate(weighted_prediction = prediction * sample_size)%>%
    group_by(site, variable) %>%
    summarise(weighted_sum = sum(weighted_prediction), 
              total_samples = sum(sample_size), .groups = "drop")%>%
    mutate(weighted_average = weighted_sum / total_samples)%>%
    select(site, variable, weighted_average)%>%
    group_by(site) %>%
    summarise(across(c(weighted_average), list(
      model_lower = ~quantile(.x, probs = 0.025),
      model_median = ~quantile(.x, probs = 0.50),
      model_upper = ~quantile(.x, probs = 0.975)
    )))
  prev_by_site <- merge(prev_by_site,fit_by_site,by="site")
  print(prev_by_site)
  names(prev_by_site) <- c("site", "observed_mean", "observed_lower", "observed_upper","Status","model_lower", "model_median", "model_upper")
  
  
  
  # plot by site setup
  site_plot <- ggplot(prev_by_site, aes(x = observed_mean, y = model_median,color=Status)) +
    geom_point() +
    geom_errorbar(aes(ymin = model_lower, ymax = model_upper), width = 0.01) +
    geom_errorbarh(aes(xmin = observed_lower, xmax = observed_upper), height = 0.01) +
    labs(x = "Observed Prevalence (Mean)", y = "Modeled Prevalence (Median)")+
    theme_minimal()+geom_abline()
  
  # Temporal analysis setup
#  temp_preds <- get_weight_month_pred(fitted_model$fit_array, survey_data)
  prev_by_month <- survey_data %>%
    group_by(month) %>%
    summarise(observed_mean = sum(positive) / sum(sample_size),
              observed_lower = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$lower,
              observed_upper = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$upper) %>%
    mutate(Status = ifelse(month %in% pred_months, "Predicted", "Fitted")) %>%
    arrange(month) %>%
    mutate(change_in_month = c(1, diff(month) > 1),
           change_in_status = c(1, diff(as.numeric(as.factor(Status))) != 0),
           Segment = cumsum(change_in_month | change_in_status)) %>%
    select(-change_in_month, -change_in_status) 
  fit_by_month<-fitted_model$data_w_probs%>%
    filter(ANC==0)%>%
    tidyr::pivot_longer(cols = starts_with("V"), names_to = "variable", values_to = "prediction") %>%
    mutate(weighted_prediction = prediction * sample_size)%>%
    group_by(month, variable) %>%
    summarise(weighted_sum = sum(weighted_prediction), 
              total_samples = sum(sample_size), .groups = "drop")%>%
    mutate(weighted_average = weighted_sum / total_samples)%>%
    select(month, variable, weighted_average)%>%
    group_by(month) %>%
    summarise(across(c(weighted_average), list(
      model_lower = ~quantile(.x, probs = 0.025),
      model_median = ~quantile(.x, probs = 0.50),
      model_upper = ~quantile(.x, probs = 0.975)
    )))
  prev_by_month <- merge(prev_by_month,fit_by_month,by="month")
  names(prev_by_month) <- c("month", "observed_mean", "observed_lower", "observed_upper","Status","Segment","model_lower", "model_median", "model_upper")
  
  #%>%
   # cbind(., rowQuantiles(temp_preds, p = c(0.025, 0.5, 0.975)))
  #names(prev_by_month) <- c("month", "observed_mean", "observed_lower", "observed_upper","Status","Segment","model_lower", "model_median", "model_upper")

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
  
  # Posterior distribution analysis
  sample_matrix <- as.matrix(fitted_model$draws)
  post_samples <- as.data.frame(sample_matrix[, !grepl("^theta", colnames(sample_matrix))])
  posterior_long <- pivot_longer(post_samples, cols = everything(), names_to = "parameter", values_to = "value")
  color_mapping <- viridis::viridis(length(unique(posterior_long$parameter)), option = "C")
  names(color_mapping) <- unique(posterior_long$parameter)
  posterior_long$color <- color_mapping[as.character(posterior_long$parameter)]
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
  
  # Composite plot
  overall_plot <- (posterior_plot | site_plot) / temporal_plot +
    plot_layout(heights = c(1, 2))
  
  return(list(
    posterior_plot = posterior_plot,
    site_plot = site_plot,
    temporal_plot = temporal_plot,
    overall_plot = overall_plot
  ))
}
