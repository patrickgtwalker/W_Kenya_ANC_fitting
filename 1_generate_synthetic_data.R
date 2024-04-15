if (!require("pacman", character.only = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(DiagrammeR, greta,dplyr,binom,matrixStats,ggplot2,bayesplot)
##install_greta_deps()
odds_to_probability <- function(odds) {
  probability <- odds / (1 + odds)
  return(probability)
}

probability_to_odds <- function(probability) {
  odds <- probability / (1 - probability)
  return(odds)
}
log_odd_ratio_site
set.seed(100)
# Data collection variables
n_month <- 50
n_sites <- 10
cMIS_samples_site_month_min <- 20
cMIS_samples_site_month_max <- 100
ANC_samples_site_month_min <- 10
ANC_samples_site_month_max <- 50
# Model parameters
cMIS_prev <- runif(n_month)  # Random uniform probabilities for each month
site_sd <- 2
log_odd_ratio_site <- rnorm(n_sites, 0, site_sd)  # Odds ratios for each site
ANC_log_odds_ratio<-rnorm(1,0,1)

# Initialize an empty dataframes for cMIS and ANC
cMIS_data <- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())
ANC_data<- data.frame(month = integer(), site = integer(), sample_size = integer(), positive = integer())
# Loop over sites and months
for (site in 1:n_sites) {
  for (month in 1:n_month) {
    # Draw a sample size for this site-month pair
    cMIS_sample_size <- sample(cMIS_samples_site_month_min:cMIS_samples_site_month_max, 1)
    ANC_sample_size <- sample(ANC_samples_site_month_min:ANC_samples_site_month_max, 1)
    
    
    # calculate site-month prevalence
    cMIS_prevalence_s_m<-odds_to_probability(probability_to_odds(cMIS_prev[month])*exp(log_odd_ratio_site[site]))
    ANC_prev_s_m<-odds_to_probability(probability_to_odds(cMIS_prevalence_s_m)*exp(ANC_log_odds_ratio))
    
    # Generate binomial variable Y
    cMIS_positive <- rbinom(1, size = cMIS_sample_size, prob = cMIS_prevalence_s_m)
    ANC_samples<-rbinom(ANC_sample_size,1,ANC_prev_s_m)
    # Add to results dataframe
    cMIS_data <- rbind(cMIS_data, c(month, site, cMIS_sample_size, cMIS_positive))
    ANC_data <- rbind(ANC_data, data.frame(month=month, site=site, sample_size=1, positive=ANC_samples))
  }
}
cMIS_prevalence_s_m
ANC_prev_s_m
rbinom(ANC_sample_size,1,ANC_prev_s_m)
# Set column names for the dataframe
colnames(cMIS_data) <- c("month", "site", "sample_size", "positive")
colnames(ANC_data) <- c("month", "site", "sample_size", "positive")

cMIS_data$ANC=0
ANC_data$ANC=1

combined_data<-rbind(ANC_data,cMIS_data)
check<-run_greta_preg_simple(combined_data,100,1000,1)
check$matrix



dim(check$matrix)
ggplot(check$prev_trends,aes(x=month,y=med))+
  geom_line()+
  geom_line(aes(y=low))+
  geom_line(aes(y=high))+
  geom_point(aes(y=cMIS_prev))



# Print the first few rows of the dataframe
cMIS_data



### data collection variables ##
n_month<-50
n_sites<-10
cMIS_samples_site_month_max<-100
cMIS_samples_site_month_max<-20


### model parameters  ##
cMIS_prev<-runif(50)
site_sd<-2
odd_ratio_site<-rnorm(10,0,site_sd)




cMIS_prev
odds_to_probability(probability_to_odds(cMIS_prev))
ANC_prev<-odds_to_probability(probability_to_odds(cMIS_prev)*exp(ANC_log_odds_ratio))
cMIS_data<-data.frame(month=1:n_month,Y=rbinom(n=n_month,size=cMIS_samples_per_month,p=cMIS_prev),N=rep(cMIS_samples_per_month,n_month),ANC=0)
ANC_samples<-unlist(
  mapply(
    function(prob,samp){
      rbinom(n=samp,p=prob,size = 1)},prob=ANC_prev,samp=ANC_samples_per_month, SIMPLIFY = FALSE
  )
)

ANC_months<-rep(seq_along(ANC_samples_per_month), times = ANC_samples_per_month)
ANC_data<-data.frame(month=ANC_months,Y=ANC_samples,N=1,ANC=1)

probability_to_odds(ANC_prev)/probability_to_odds(cMIS_prev)
exp(ANC_log_odds_ratio)

check2<-ANC_data%>%
  group_by(month)%>%
  summarise(prev=sum(Y)/sum(N))

plot(check2$prev,ANC_prev)

ANC_log_odds_ratio
mcmc_intervals(check$draws)


as.matrix(calculate(1/(1+exp(-theta_month)),values=draws))
       check$draws
       plot(check$prev_trends)
       
  check$prev_trends
combined_data
rbinom(n=rep(ANC_samples_per_month,length(ANC_prev)),p=ANC_prev,size = 1)
#<-data.frame(month=1:n_month,Y=rbinom(samples_per_month,n_month,cMIS_prev),N=rep(samples_per_month,n_month))
rbinom(n=1,p=c(0.2,0.5),size=c(10,20))
rbinom(n=10,p=0.2,size=1)


ANC_data<-data.frame()

