source("setup.R")
source("functions.R")

set.seed(100)
# Data collection variables
n_month <- 50
n_sites <- 10
cMIS_samples_site_month_min <- 20
cMIS_samples_site_month_max <- 100
ANC_samples_site_month_min <- 10
ANC_samples_site_month_max <- 50
pred_months<-c(4:28,32:50)


source("data_generation_simple.R")
simple_data<-generate_data_simple(n_month,
                               n_sites,
                               cMIS_samples_site_month_min,
                               cMIS_samples_site_month_max,
                               ANC_samples_site_month_min,
                               ANC_samples_site_month_max,
                               pred_months
                               
)
source("models.R")
simple_fitting<-run_greta_preg_simple(simple_data$data_to_model,100,200,1)

source("plot_functions.R")
plots<-get_plots(simple_fitting,simple_data$simulated_survey_data,simple_data$param_df,pred_months)
plots
