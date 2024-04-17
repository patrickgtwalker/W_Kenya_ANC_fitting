
get_plots <- function(fitted_model, survey_data, param_df, pred_months=F, pred_sites=F) {
 
  cols<-viridis(2)
  
  # Extract predictions and observed statistics
  site_preds <- get_weight_site_pred(fitted_model$fit_array, survey_data)
  prev_by_site <- survey_data %>%
    group_by(site) %>%
    summarise(observed_mean = sum(positive) / sum(sample_size),
              observed_lower = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$lower,
              observed_upper = binom.confint(sum(positive), sum(sample_size), methods = "wilson")$upper)%>%
    mutate(Status = ifelse(site %in% pred_sites, "Predicted", "Fitted"))
  prev_by_site <- cbind(prev_by_site, rowQuantiles(site_preds, p = c(0.025, 0.5, 0.975)))
  names(prev_by_site) <- c("site", "observed_mean", "observed_lower", "observed_upper","Status","model_lower", "model_median", "model_upper")
  
  # Spatial plot setup
  spatial_plot <- ggplot(prev_by_site, aes(x = observed_mean, y = model_median,color=Status)) +
    geom_point() +
    geom_errorbar(aes(ymin = model_lower, ymax = model_upper), width = 0.01) +
    geom_errorbarh(aes(xmin = observed_lower, xmax = observed_upper), height = 0.01) +
    labs(x = "Observed Prevalence (Mean)", y = "Modeled Prevalence (Median)")+
    theme_minimal()+geom_abline()
  spatial_plot
  # Temporal analysis setup
  temp_preds <- get_weight_month_pred(fitted_model$fit_array, survey_data)
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
    select(-change_in_month, -change_in_status) %>%
    cbind(., rowQuantiles(temp_preds, p = c(0.025, 0.5, 0.975)))
  names(prev_by_month) <- c("month", "observed_mean", "observed_lower", "observed_upper","Status","Segment","model_lower", "model_median", "model_upper")

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
  temporal_plot
  # Posterior distribution analysis
  sample_matrix <- as.matrix(fitted_model$draws)
  post_samples <- as.data.frame(sample_matrix[, !grepl("^theta", colnames(sample_matrix))])
  posterior_long <- pivot_longer(post_samples, cols = everything(), names_to = "parameter", values_to = "value")
  color_mapping <- viridis::viridis(length(unique(posterior_long$parameter)), option = "C")
  names(color_mapping) <- unique(posterior_long$parameter)
  posterior_long$color <- color_mapping[as.character(posterior_long$parameter)]
  simulated_long <- pivot_longer(param_df, cols = everything(), names_to = "parameter", values_to = "value")
  simulated_long$color <- color_mapping[as.character(simulated_long$parameter)]
  print(color_mapping)
  print(simulated_long)
  posterior_plot <- ggplot(posterior_long, aes(x = value, y = parameter, fill = color)) +
    geom_density_ridges(rel_min_height = 0.01, alpha = 0.25) +
    scale_fill_identity() +
    geom_vline(data = simulated_long, aes(xintercept = value, color = color), linetype = "dashed", size = 1) +
    scale_color_identity() +
    theme_ridges() +
    labs(title = "Posterior Distributions with Simulated Values", x = "Parameter Values", y = "Parameter") +
    theme(legend.position = "none")
  
  # Composite plot
  overall_plot <- (posterior_plot | spatial_plot) / temporal_plot +
    plot_layout(heights = c(1, 2))
  
  return(list(
    posterior_plot = posterior_plot,
    spatial_plot = spatial_plot,
    temporal_plot = temporal_plot,
    overall_plot = overall_plot
  ))
}
