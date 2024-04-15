
run_greta_preg_simple<-function(data,warmup,n_samples,thin){
  months=max(data$month)
  sites=max(data$site)
  theta_month<-normal(0,100,months)
  tau_site<-gamma(0.01,0.01)
  theta_site<-normal(0,tau_site,sites)
  preg_par<-normal(0,100)
  p<-ilogit(preg_par*data$ANC+theta_month[data$month]+theta_site[data$site])
  distribution(data$positive)=binomial(data$sample_size,p)
  m<-model(preg_par,theta_month,theta_site,tau_site)
  draws <- mcmc(m, n_samples = n_samples, warmup = warmup,thin=thin)
  fit_array<-array(0,dim = c(n_samples/thin*4,sites,months))
  for(i in 1:sites){
    cmis_log_odds_month_site<-theta_site[i]+theta_month
    fit_array[,i,]<-as.matrix(calculate(odds_to_probability(exp(cmis_log_odds_month_site)),values=draws))
  }
  matrix<-as.matrix(calculate(odds_to_probability(exp(theta_month)),values=draws))
  prev_trends<-as.data.frame(colQuantiles(matrix,probs=c(0.5,0.025,0.975)))
  prev_trends$month=1:months
  names(prev_trends)<-c("med","low","high","month")
  return(
    list(
      data=data,
      draws=draws,
      prev_trends=prev_trends,
      fit_array=fit_array
    )
  )
}
