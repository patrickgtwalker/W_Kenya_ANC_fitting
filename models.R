deviance <- function(prob, y, n){
  return(-2 *sum(dbinom(y, size = n, prob = prob, log = TRUE)))
}

DIC_calc<-function(data,probs,REs_site_dev,REs_month_dev,p_temp,site_ave,tau_site_ave,month_ave,tau_month_ave){
  logit_dev <- unlist(lapply(X = probs,
                             FUN = deviance,
                             y = data$positive,
                             n = data$sample_size))
  p_ave <- as.data.frame(do.call(rbind, calculate(p_temp)))
  Re_site_hat<- -2*sum(log(1 / (tau_site_ave * sqrt(2 * pi)) * exp(-((site_ave) / tau_site_ave) ^ 2/2)))
  Re_month_hat<- -2*sum(log(1 / (tau_month_ave * sqrt(2 * pi)) * exp(-((month_ave) / tau_month_ave) ^ 2/2)))
  logit_hat <- deviance(prob=p_ave$V1,
                        y = data$positive,
                        n = data$sample_size)
  mean_dev<-quantile(logit_dev+REs_site_dev+REs_month_dev,p=0.5)
  d_hat<-logit_hat+Re_site_hat+Re_month_hat
  #p_d is the effective number of parameters
  p_d <- mean_dev - d_hat
  DIC <- p_d + mean_dev
  return(list(
    DIC=DIC,
    p_d=p_d
  ))
}

fit_m1<-function(data,warmup,n_samples,thin){
  months=max(data$month)
  sites=max(data$site)
  tau_month<-gamma(0.01,0.01)
  tau_site<-gamma(0.01,0.01)
  theta_month<-normal(0,tau_month,months)
  theta_site<-normal(0,tau_site,sites)
  base_par<-normal(0,100)
  preg_par<-normal(0,100)
  p<-ilogit(base_par+preg_par*data$ANC+theta_month[data$month]+theta_site[data$site])
  distribution(data$positive)=binomial(data$sample_size,p)
  m<-model(base_par,preg_par,theta_month,tau_month,theta_site,tau_site)
  draws <- mcmc(m, n_samples = n_samples, warmup = warmup,thin=thin)
  #fit_array<-array(0,dim = c(n_samples/thin*4,sites,months))
  #for(i in 1:sites){
  #  cmis_log_odds_month_site<-calculate(base_par+theta_site[i]+theta_month,values=draws)
   # fit_array[,i,]<-odds_to_probability(exp(as.matrix(cmis_log_odds_month_site)))
  #}
  #matrix<-odds_to_probability(exp(as.matrix(calculate(base_par+theta_month,values=draws))))
  #prev_trends<-as.data.frame(colQuantiles(matrix,probs=c(0.5,0.025,0.975)))
  #prev_trends$month=1:months
  #names(prev_trends)<-c("med","low","high","month")
  probs <- as.data.frame(t(do.call(rbind, calculate(p, values=draws))))
  data_w_probs<-cbind(data,probs)
  lik_site <- (1 / (tau_site * sqrt(2 * pi))) * exp(-((theta_site) / tau_site) ^ 2/2)
  lik_month <- (1 / (tau_month * sqrt(2 * pi))) * exp(-((theta_month) / tau_month) ^ 2/2)
  REs_site_dev<-unlist(calculate(-2*sum(log(lik_site)),values=draws))
  REs_month_dev<-unlist(calculate(-2*sum(log(lik_month)),values=draws))
  
  median_coefs <- as.list(colQuantiles(do.call(rbind,draws),p=0.5))
  base_ave <- unlist(median_coefs[grep('base_par',names(median_coefs))])
  preg_ave <- unlist(median_coefs[grep('preg_par',names(median_coefs))])
  month_ave <- unlist(median_coefs[grep('theta_month',names(median_coefs))])
  site_ave <- unlist(median_coefs[grep('theta_site',names(median_coefs))])
  tau_site_ave <- unlist(median_coefs[grep('tau_site',names(median_coefs))])
  tau_month_ave <- unlist(median_coefs[grep('tau_month',names(median_coefs))])
  p_temp <- ilogit(base_ave+preg_ave*data$ANC+month_ave[data$month]+site_ave[data$site])
  #Use calculate to turn the greta array into a dataframe
  DIC_calc=DIC_calc(data,probs,REs_site_dev,REs_month_dev,p_temp,site_ave,tau_site_ave,month_ave,tau_month_ave)
  
  
  return(
    list(
      draws=draws,
   #   prev_trends=prev_trends,
    #  fit_array=fit_array,
      data_w_probs=data_w_probs,
      p_d=DIC_calc$p_d,
      DIC=DIC_calc$DIC
    )
  )
}

fit_m5<-function(data,warmup,n_samples,thin){
  months=max(data$month)
  sites=max(data$site)
  tau_month<-gamma(0.01,0.01)
  tau_site<-gamma(0.01,0.01)
  base_par<-normal(0,100)
  theta_month<-normal(0,tau_month,months)
  theta_site<-normal(0,tau_site,sites)
  log_odd_ratio_intercept_GC1<-normal(0,100)
  log_odd_ratio_intercept_GC2<-normal(0,100)
  log_odd_ratio_intercept_GC3<-normal(0,100)
  log_odd_ratio_gradient_GC1<-normal(0,100)
  log_odd_ratio_gradient_GC2<-normal(0,100)
  log_odd_ratio_gradient_GC3<-normal(0,100)
  p<-ilogit(base_par+theta_site[data$site]+theta_month[data$month] +
              (log_odd_ratio_intercept_GC1+theta_site[data$site]*log_odd_ratio_gradient_GC1)*data$gravcat1+
              (log_odd_ratio_intercept_GC2+theta_site[data$site]*log_odd_ratio_gradient_GC2)*data$gravcat2+
              (log_odd_ratio_intercept_GC3+theta_site[data$site]*log_odd_ratio_gradient_GC3)*data$gravcat3
  )
  distribution(data$positive)=binomial(data$sample_size,p)
  m<-model(base_par,theta_month,tau_month,theta_site,tau_site,log_odd_ratio_gradient_GC1,log_odd_ratio_gradient_GC2,log_odd_ratio_gradient_GC3,log_odd_ratio_intercept_GC1,log_odd_ratio_intercept_GC2,log_odd_ratio_intercept_GC3)
  draws <- mcmc(m, n_samples = n_samples, warmup = warmup,thin=thin)
  fit_array<-array(0,dim = c(n_samples/thin*4,sites,months))
  for(i in 1:sites){
    cmis_log_odds_month_site<-calculate(base_par+theta_site[i]+theta_month,values=draws)
    fit_array[,i,]<-odds_to_probability(exp(as.matrix(cmis_log_odds_month_site)))
  }
  matrix<-odds_to_probability(exp(as.matrix(calculate(base_par+theta_month,values=draws))))
  prev_trends<-as.data.frame(colQuantiles(matrix,probs=c(0.5,0.025,0.975)))
  prev_trends$month=1:months
  names(prev_trends)<-c("med","low","high","month")
  probs <- as.data.frame(t(do.call(rbind, calculate(p, values=draws))))
  lik_site <- (1 / (tau_site * sqrt(2 * pi))) * exp(-((theta_site) / tau_site) ^ 2/2)
  lik_month <- (1 / (tau_month * sqrt(2 * pi))) * exp(-((theta_month) / tau_month) ^ 2/2)
  REs_site_dev<-unlist(calculate(-2*sum(log(lik_site)),values=draws))
  REs_month_dev<-unlist(calculate(-2*sum(log(lik_month)),values=draws))
  
  median_coefs <- as.list(colQuantiles(do.call(rbind,draws),p=0.5))
  base_ave <- unlist(median_coefs[grep('base_par',names(median_coefs))])
  log_odd_ratio_intercept_GC1_ave<-unlist(median_coefs[grep('log_odd_ratio_intercept_GC1',names(median_coefs))])
  log_odd_ratio_intercept_GC2_ave<-unlist(median_coefs[grep('log_odd_ratio_intercept_GC2',names(median_coefs))])
  log_odd_ratio_intercept_GC3_ave<-unlist(median_coefs[grep('log_odd_ratio_intercept_GC3',names(median_coefs))])
  log_odd_ratio_gradient_GC1_ave<-unlist(median_coefs[grep('log_odd_ratio_gradient_GC1',names(median_coefs))])
  log_odd_ratio_gradient_GC2_ave<-unlist(median_coefs[grep('log_odd_ratio_gradient_GC2',names(median_coefs))])
  log_odd_ratio_gradient_GC3_ave<-unlist(median_coefs[grep('log_odd_ratio_gradient_GC3',names(median_coefs))])
  month_ave <- unlist(median_coefs[grep('theta_month',names(median_coefs))])
  site_ave <- unlist(median_coefs[grep('theta_site',names(median_coefs))])
  tau_site_ave <- unlist(median_coefs[grep('tau_site',names(median_coefs))])
  tau_month_ave <- unlist(median_coefs[grep('tau_month',names(median_coefs))])
  p_temp <-  ilogit(base_ave+site_ave[data$site]+month_ave[data$month] +
           (log_odd_ratio_intercept_GC1_ave+site_ave[data$site]*log_odd_ratio_gradient_GC1_ave)*data$gravcat1+
           (log_odd_ratio_intercept_GC2_ave+site_ave[data$site]*log_odd_ratio_gradient_GC2_ave)*data$gravcat2+
           (log_odd_ratio_intercept_GC3_ave+site_ave[data$site]*log_odd_ratio_gradient_GC3_ave)*data$gravcat3
  )
  
  #Use calculate to turn the greta array into a dataframe
  DIC_calc=DIC_calc(data,probs,REs_site_dev,REs_month_dev,p_temp,site_ave,tau_site_ave,month_ave,tau_month_ave)
  
  
  return(
    list(
      data=data,
      draws=draws,
      prev_trends=prev_trends,
      fit_array=fit_array,
      probs=probs,
      p_d=as.numeric(DIC_calc$p_d),
      DIC=as.numeric(DIC_calc$DIC)
    )
  )
}

predict_months_m1<-function(data,warmup,n_samples,thin){
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


compare_RE_site<-function(data,warmup,n_samples,thin){
  
}

predict_sites_m1<-function(data,warmup,n_samples,thin){
  months=max(data$month)
  sites=max(data$site)
  tau_month<-gamma(0.01,0.01)
  theta_month<-normal(0,tau_month,months)
  theta_site<-normal(0,100,sites)
  preg_par<-normal(0,100)
  p<-ilogit(preg_par*data$ANC+theta_month[data$month]+theta_site[data$site])
  distribution(data$positive)=binomial(data$sample_size,p)
  m<-model(preg_par,theta_month,theta_site,tau_month)
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
