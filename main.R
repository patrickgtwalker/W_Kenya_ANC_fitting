source("setup.R")
source("functions.R")
source("data_generation.R")
source("models.R")
source("plot_functions.R")
complex_data<-generate_data_complex()
complex_fitting<-run_greta_preg_complex(complex_data$data_to_model,100,200,1)
`complex_fitting$p_d

simple_fitting<-run_greta_preg_simple(simple_data$data_to_model,100,200,1)
set.seed(100)
get_plots(predict_spatial_data,site_data_to_predict$simulated_survey_data,site_data_to_predict$param_df,pred_sites=site_data_to_predict$pred_sites)
simple_fitting$DIC

simple_data<-generate_data_simple()


simple_data$data_to_model
simple_data$
simple_fitting<-run_greta_preg_simple(simple_data$data_to_model,100,200,1)
complex_fitting<-run_greta_preg_complex(complex_data$data_to_model,100,200,1)
simple_fitting_complex_data<-run_greta_preg_simple(complex_data$data_to_model,100,200,1)

simple_plots<-get_plots(simple_fitting,simple_data$simulated_survey_data,simple_data$param_df)
complex_plots<-get_plots(complex_fitting,complex_data$simulated_survey_data,complex_data$param_df)
simple_fitting_complex_data_plot<-get_plots(simple_fitting_complex_data,complex_data$simulated_survey_data,complex_data$param_df)




temporal_data_to_predict<-generate_data_simple_NNE()
predict_temp_data<-run_greta_predict_month(temporal_data_to_predict$data_to_model,100,200,1)

site_data_to_predict<-generate_data_simple_FE_site()
check<-site_data_to_predict$data_to_model
predict_spatial_data<-run_greta_predict_site(site_data_to_predict$data_to_model,1000,2000,1)
check<-site_data_to_predict$data_to_model
checksite_data_to_predict$pred_sites

check<-site_data_to_predict$data_to_model


mcmc_intervals(predict_spatial_data$draws)

simple_fitting_complex_data_plot
complex_plots
mcmc_intervals(simple_fitting$draws)
complex_plots
plots$temporal_plot
simple_plots$spatial_plot
dim(complex_fitting$fit_array)
mcmc_intervals(complex_fitting$draws)
rbinom(4, size = 1000, prob = c(0.1,0.2,0.3,0.4))
complex_fitting$simulated_survey_data
simple_fitting$
  
  prev_site<-complex_data$data_to_model%>%group_by(site,ANC)%>%
  summarise(prev=sum(positive)/sum(sample_size))%>%
  pivot_wider(names_from = ANC, values_from = prev, names_prefix = "ANC_")

ggplot(prev_site, aes(x = ANC_0, y = ANC_1)) +
  geom_point() +  # Add points
  geom_text(aes(label = site), vjust = -1) +  # Add site labels
  labs(x = "Prevalence when ANC=0", y = "Prevalence when ANC=1", title = "Site-specific Prevalence Comparison") +
  theme_minimal() +  # Use a minimal theme
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") 
simple_data$param_df