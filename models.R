deviance <- function(prob, y, n){
  return(-2 *sum(dbinom(y, size = n, prob = prob, log = TRUE)))
}

DIC_calc<-function(data,probs,REs_site_dev,REs_month_dev,p_temp,site_ave,tau_site_ave,month_ave,tau_month_ave){
  print("here1")
  logit_dev <- unlist(lapply(X = probs,
                             FUN = deviance,
                             y = data$positive,
                             n = data$sample_size))
  print("here2")
  p_ave <- as.data.frame(do.call(rbind, calculate(p_temp)))
  Re_site_hat<- -2*sum(log(1 / (tau_site_ave * sqrt(2 * pi)) * exp(-((site_ave) / tau_site_ave) ^ 2/2)))
  Re_month_hat<- -2*sum(log(1 / (tau_month_ave * sqrt(2 * pi)) * exp(-((month_ave) / tau_month_ave) ^ 2/2)))
  print("here3")
  logit_hat <- deviance(prob=p_ave$V1,
                        y = data$positive,
                        n = data$sample_size)
  mean_dev<-quantile(logit_dev+REs_site_dev+REs_month_dev,p=0.5)
  d_hat<-logit_hat+Re_site_hat+Re_month_hat
  #p_d is the effective number of parameters
  print("here4")
  p_d <- mean_dev - d_hat
  DIC <- p_d + mean_dev
  return(list(
    DIC=DIC,
    p_d=p_d
  ))
}

run_greta_preg_simple<-function(data,warmup,n_samples,thin){
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
  fit_array<-array(0,dim = c(n_samples/thin*4,sites,months))
  for(i in 1:sites){
    cmis_log_odds_month_site<-base_par+theta_site[i]+theta_month
    fit_array[,i,]<-as.matrix(calculate(odds_to_probability(exp(cmis_log_odds_month_site)),values=draws))
  }
  matrix<-as.matrix(calculate(odds_to_probability(exp(theta_month)),values=draws))
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
      data=data,
      draws=draws,
      prev_trends=prev_trends,
      fit_array=fit_array,
      p_d=DIC_calc$p_d,
      DIC=DIC_calc$DIC
    )
  )
}

run_greta_preg_complex<-function(data,warmup,n_samples,thin){
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
    cmis_log_odds_month_site<-theta_site[i]+theta_month
    fit_array[,i,]<-as.matrix(calculate(odds_to_probability(exp(cmis_log_odds_month_site)),values=draws))
  }
  matrix<-as.matrix(calculate(odds_to_probability(exp(theta_month)),values=draws))
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
      p_d=DIC_calc$p_d,
      DIC=DIC_calc$DIC
    )
  )
}

run_greta_predict_month<-function(data,warmup,n_samples,thin){
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


run_greta_predict_site<-function(data,warmup,n_samples,thin){
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
run_greta_preg_casdc<-function(data,warmup,n_samples,thin){
  data$month_code<-round(data$month_code)
  data<-data%>%
    mutate(primi=ifelse(grav==1,1,0))%>%
    mutate(primi=replace(primi,is.na(primi),0))%>%
    mutate(multi=ifelse(grav%in%c(2,3),1,0))%>%
    mutate(grand=ifelse(grav>=4,1,0))%>%
    mutate(grand=replace(grand,is.na(grand),0))
  site_col<-grep(site_str,names(data))
  sites<-max(data[,site_col])
  months<-round(max(data[,grep("month_code",names(data))]))
  
  tau_month<-gamma(0.01,0.01)
  tau_site<-gamma(0.01,0.01)
  theta_month<-normal(0,tau_month,months)
  theta_site<-normal(0,tau_site,sites)
  base_par<-normal(0,100)
  log_odd_ratio_intercept_GC1<-normal(0,100)
  log_odd_ratio_intercept_GC2<-normal(0,100)
  log_odd_ratio_intercept_GC3<-normal(0,100)
  log_odd_ratio_gradient_GC1<-normal(0,100)
  log_odd_ratio_gradient_GC2<-normal(0,100)
  log_odd_ratio_gradient_GC3<-normal(0,100)
  p<-ilogit(base_par+trans_par*theta_site[as.numeric(adjust_data$site)]*(adjust_data$primi+adjust_data$multi+adjust_data$grand)+primi_par*adjust_data$primi+multi_par*adjust_data$multi+grand_par*adjust_data$grand+theta_month[as.numeric(adjust_data$month_code)]+theta_site[as.numeric(adjust_data$site)])
  distribution(adjust_data$Y)=binomial(adjust_data$N,p)
  m<-model(base_par,trans_par,primi_par,multi_par,grand_par,theta_month,tau_month,theta_site,tau_site)
  draws <- mcmc(m, n_samples = n_samples, warmup = warmup,thin=thin)
  probs <- as.data.frame(t(do.call(rbind, calculate(p, values=draws))))
  logit_dev <- unlist(lapply(X = probs,
                             FUN = deviance,
                             y = adjust_data$Y,
                             n = adjust_data$N))
  
  lik_site <- (1 / (tau_site * sqrt(2 * pi))) * exp(-((theta_site) / tau_site) ^ 2/2)
  lik_month <- (1 / (tau_month * sqrt(2 * pi))) * exp(-((theta_month) / tau_month) ^ 2/2)
  REs_site_dev<-unlist(calculate(-2*sum(log(lik_site)),values=draws))
  REs_month_dev<-unlist(calculate(-2*sum(log(lik_month)),values=draws))
  # Calculate mean of each coefficient
  mean_coefs <- as.list(colQuantiles(do.call(rbind,draws),p=0.5))
  
  #Calculate p based on the mean coefficients
  #The first lines are repackaging the mean coefficients back into a structure that works with the greta calculate function.
  base_ave <- unlist(mean_coefs[grep('base_par',names(mean_coefs))])
  trans_par_ave <- unlist(mean_coefs[grep('trans_par',names(mean_coefs))])
  primi_ave <- unlist(mean_coefs[grep('primi_par',names(mean_coefs))])
  multi_ave <- unlist(mean_coefs[grep('multi_par',names(mean_coefs))])
  grand_ave <- unlist(mean_coefs[grep('grand_par',names(mean_coefs))])
  month_ave <- unlist(mean_coefs[grep('theta_month',names(mean_coefs))])
  site_ave <- unlist(mean_coefs[grep('theta_site',names(mean_coefs))])
  tau_site_ave <- unlist(mean_coefs[grep('tau_site',names(mean_coefs))])
  tau_month_ave <- unlist(mean_coefs[grep('tau_month',names(mean_coefs))])
  
  #p_temp is the greta array for the regression model but using the mean coefficients
  p_temp <- ilogit(base_ave+trans_par_ave*site_ave[as.numeric(adjust_data$site)]*(adjust_data$primi+adjust_data$multi+adjust_data$grand)+primi_ave*adjust_data$primi+multi_ave*adjust_data$multi+grand_ave*adjust_data$grand+month_ave[adjust_data$month_code]+site_ave[adjust_data$site])
  #Use calculate to turn the greta array into a dataframe
  p_ave <- as.data.frame(do.call(rbind, calculate(p_temp)))
  Re_site_hat<- -2*sum(log(1 / (tau_site_ave * sqrt(2 * pi)) * exp(-((site_ave) / tau_site_ave) ^ 2/2)))
  Re_month_hat<- -2*sum(log(1 / (tau_month_ave * sqrt(2 * pi)) * exp(-((month_ave) / tau_month_ave) ^ 2/2)))
  logit_hat <- deviance(prob=p_ave$V1,
                        y = adjust_data$Y,
                        n = adjust_data$N)
  mean_dev<-quantile(logit_dev+REs_site_dev+REs_month_dev,0.5)
  d_hat<-logit_hat+Re_site_hat+Re_month_hat
  #p_d is the effective number of parameters
  p_d <- mean_dev - d_hat
  DIC <- p_d + mean_dev
  
  #Log-likelihood
  #ll <- dev/-2
  
  return(list(draws=draws,
              mean_dev=mean_dev,
              mean_logit=mean(logit_dev),
              REs_site_dev=mean(REs_site_dev),
              REs_month_dev=mean(REs_month_dev),
              mean_RE_site=mean(REs_site_dev),
              mean_RE_month=mean(REs_month_dev),
              d_hat=d_hat,
              logit_hat=logit_hat,
              Re_site_hat=Re_site_hat,
              Re_month_hat=Re_month_hat,
              p_d=p_d,
              DIC=DIC
  ))
}
