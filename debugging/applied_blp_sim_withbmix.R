## ----------------------------- Build Functions -----------------------------------##
# TLP
wrapper_TLP <- function(component_forecasts, component_pits, weights, eval = TRUE, pit = FALSE){
  if(eval){
    weights<-exp(weights)
    weights <- weights/sum(weights)
  }
  ll <- log(component_forecasts%*% weights)
  dens <- component_forecasts%*% weights
  PIT <- component_pits %*% weights
  ls <- ifelse(eval,sum(ll),mean(ll))
  if(pit){
    return(list(ls,PIT,dens))
  } else{
    return(ls)
  }
}

# BLP
wrapper_BLP <- function(component_forecasts, component_pits, pars, eval = TRUE, pit = FALSE){
  if(eval){
    ab     <- exp(pars[1:2])
    weights <- exp(pars[-c(1,2)])
    weights <- weights/sum(weights)
  }else{
    ab     <- pars[1:2]
    weights <- pars[-c(1,2)]
  }
  g      <- component_forecasts %*% weights
  g_beta <- dbeta(component_pits %*% weights, shape1 = ab[1], shape2 = ab[2])
  dens <- g*g_beta
  ll     <- log(g) + log(g_beta)
  PIT <- pbeta(component_pits %*% weights, shape1 = ab[1], shape2 = ab[2])
  ls <- ifelse(eval,sum(ll),mean(ll))
  if(pit){
    return(list(ls,PIT,dens))
  } else{
    return(ls)
  }
} 
# EW
wrapper_EW <- function(component_forecasts, component_pits, eval = TRUE, pit = FALSE){
  weights <- rep(1/ncol(component_forecasts),ncol(component_forecasts))
  ll <- log(component_forecasts %*% weights)
  dens <- component_forecasts %*% weights
  PIT <- component_pits %*% weights
  ls <- ifelse(eval,sum(ll),mean(ll))
  if(pit){
    return(list(ls,PIT,dens))
  } else{
    return(ls)
  }
}

# EW-BLP
wrapper_ewBLP <- function(component_forecasts, component_pits, pars, eval = TRUE, pit = FALSE){
  if(eval){
    ab     <- exp(pars)
    weights <- exp(rep(1/ncol(component_forecasts),ncol(component_forecasts)))
    weights <- weights/sum(weights)
  }else{
    ab     <- pars
    weights <- rep(1/ncol(component_forecasts),ncol(component_forecasts))
  }
  g      <- component_forecasts %*% weights
  g_beta <- dbeta(component_pits %*% weights,shape1 = ab[1], shape2 = ab[2])
  dens <- g*g_beta
  ll     <- log(g) + log(g_beta)
  PIT <- pbeta(component_pits %*% weights, shape1 = ab[1], shape2 = ab[2])
  ls <- ifelse(eval,sum(ll),mean(ll))
  if(pit){
    return(list(ls,PIT,dens))
  } else{
    return(ls)
  }
}
##------------------------------Result Functions -----------------------------------##

get_component_res <- function(component_forecasts,component_pits,train_id, test_id,PIT=TRUE){
  train_logscore <- log(component_forecasts[train_id,i])
  test_logscore <- log(component_forecasts[test_id,i])
  train_pits <- component_pits[train_id,i]
  test_pits <- component_pits[test_id,i]
  if(PIT){
    results <- list(train_logscore,test_logscore,train_pits,test_pits)
    names(results) <- c("train_logscore","test_logscore","train_pits","test_pits")
    } 
  else{
    results <- list(train_logscore,test_logscore)
    names(results) <- c("train_logscore","test_logscore")
    }
  return(results)
}

# make frequentist ensembles
make_freq_ensemble <- function(ensemble_name,component_forecasts,component_pits, initial_vals, train_id, test_id){
  results <- list()
  # TLP 
  if(ensemble_name=="TLP"){
    # get "params"
    weights_TLP <- log(initial_vals)
    est <- optim(
      par     = weights_TLP,
      fn      = wrapper_TLP,
      component_forecasts = component_forecasts[train_id, ],
      component_pits = component_pits[train_id, ],
      method  = "Nelder-Mead",
      control = list(fnscale = -1, trace = TRUE, maxit = 500)
    )
    weights_TLP <- exp(est$par)/sum(exp(est$par))
    results[["params"]] <- weights_TLP
    train_res <- wrapper_TLP(component_forecasts[train_id, ], component_pits[train_id, ], results$params, eval = FALSE, pit = TRUE)
    test_res <- wrapper_TLP(component_forecasts[test_id, ], component_pits[test_id, ], results$params, eval = FALSE, pit = TRUE)
  }
  
  # BLP 
  if(ensemble_name=="BLP"){
    # get "params"
    weights_BLP <- log(initial_vals)
    ab          <- log(c(1,1))
    pars_BLP    <- c(ab, weights_BLP)
    est <- optim(
      par     = pars_BLP,
      fn      = wrapper_BLP,
      component_forecasts = component_forecasts[train_id, ],
      component_pits = component_pits[train_id, ],
      method  = "Nelder-Mead",
      control = list(fnscale = -1, trace = TRUE, maxit = 500)
      )
    ab          <- exp(est$par[1:2])
    weights_BLP <- exp(est$par[-c(1,2)])/sum(exp(est$par[-c(1,2)]))
    pars_BLP    <- c(ab, weights_BLP)
    results[["params"]] <- pars_BLP
    train_res <- wrapper_BLP(component_forecasts[train_id, ], component_pits[train_id, ], results$params, eval = FALSE, pit = TRUE)
    test_res <- wrapper_BLP(component_forecasts[test_id, ], component_pits[test_id, ], results$params, eval = FALSE, pit = TRUE)
  }
  
  # EW
  if(ensemble_name=="EW"){
    results[["params"]] <- rep(1/ncol(component_forecasts),ncol(component_forecasts))
    train_res <- wrapper_EW(component_forecasts[train_id, ], component_pits[train_id, ], eval = FALSE, pit = TRUE)
    test_res <- wrapper_EW(component_forecasts[test_id, ], component_pits[test_id, ], eval = FALSE, pit = TRUE)
  }
  
  # EW-BLP
  if(ensemble_name=="EW-BLP"){
    pars_ewBLP  <- log(c(1,1))
    est <- optim(
      par     = pars_ewBLP,
      fn      = wrapper_ewBLP,
      component_forecasts = component_forecasts[train_id, ],
      component_pits = component_pits[train_id, ],
      method  = "Nelder-Mead",
      control = list(fnscale = -1, trace = TRUE, maxit = 500)
      )
    pars_ewBLP    <- exp(est$par)
    results[["params"]] <- pars_ewBLP
    train_res <- wrapper_ewBLP(component_forecasts[train_id, ], component_pits[train_id, ], results$params, eval = FALSE, pit = TRUE)
    test_res <- wrapper_ewBLP(component_forecasts[test_id, ], component_pits[test_id, ], results$params, eval = FALSE, pit = TRUE)
  }
  results <- c(results,train_res,test_res)
  names(results) <- c("params","train_ls","train_pits","train_dens","test_ls","test_pits","test_dens")
  return(results)
}

# make mixture ensembles
make_mix_ensemble <- function(ensemble_name,K,component_forecasts,component_pits, initial_vals,hyperpars,
                              train_id, test_id){
  M <- ncol(component_forecasts)
  if((!grepl("EW",ensemble_name,fixed = TRUE)) && (!grepl("BMC1",ensemble_name,fixed = TRUE))){
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> K;
    int<lower=0> M;
    vector<lower=0>[K] alpha_w;
    vector<lower=0>[M] alpha_omega;
    matrix[n,M] comp_forecasts;
    matrix[n,M] comp_pits;
    }
    
    transformed data {
    matrix[n,M] log_f = log(comp_forecasts);
    }
    
    parameters {
    vector<lower=0, upper=1>[K] mu;
    vector<lower=0>[K] nu;
    simplex[M] omega_k[K];
    simplex[K] w;
    }
    
    model {
    matrix[n, K] H;
    vector[K] log_weighted_mixt_density;
    vector[M] inner_log_mixt_density;
    vector[K] log_w = log(w);
    vector[M] log_omega_k[K];
    // priors
    for (k in 1:K) {
    mu[k] ~ beta(2, 2);
    nu[k] ~ gamma(0.1, 0.1);
    omega_k[k] ~ dirichlet(alpha_omega);
    log_omega_k[k] = log(omega_k[k]);
    H[, k] = comp_pits*omega_k[k];
    }
    w ~ dirichlet(alpha_w);
    for (i in 1:n) {
    log_weighted_mixt_density = log_w;
    for (k in 1:K) {
    inner_log_mixt_density = log_omega_k[k];
    for (m in 1:M) {
    inner_log_mixt_density[m] += log_f[i, m];
    }
    log_weighted_mixt_density[k] += log_sum_exp(inner_log_mixt_density);
    log_weighted_mixt_density[k] += beta_proportion_lpdf(H[i, k] | mu[k], nu[k]);
    }
    target += log_sum_exp(log_weighted_mixt_density);
    }
    }
    "
    input_data <- list(n=sum(ind_est), K=K, M=M, comp_forecasts=component_forecasts[train_id,], comp_pits=component_pits[train_id,],
                        alpha_w=hyperpars$alpha_w,alpha_omega=hyperpars$alpha_omega)
  }
  if(grepl("EW",ensemble_name,fixed = TRUE) && (!grepl("BMC1",ensemble_name,fixed = TRUE))){
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> K;
    int<lower=0> M;
    vector<lower=0>[K] alpha_w;
    vector<lower=0, upper=1>[M] omega;
    matrix[n,M] comp_forecasts;
    matrix[n,M] comp_pits;
  }
    
    transformed data {
    matrix[n,M] log_f = log(comp_forecasts);
    vector[M] log_omega = log(omega);
    vector[n] H = comp_pits*omega;
    }
    
    parameters {
    vector<lower=0, upper=1>[K] mu;
    vector<lower=0>[K] nu;
    simplex[K] w;
    }
    
    model {
    vector[K] log_weighted_mixt_density;
    vector[M] inner_log_mixt_density;
    vector[K] log_w = log(w);
    // priors
    for (k in 1:K) {
    mu[k] ~ beta(2, 2);
    nu[k] ~ gamma(0.1, 0.1);
    }
    w ~ dirichlet(alpha_w);
    for (i in 1:n) {
    log_weighted_mixt_density = log_w;
    for (k in 1:K) {
    inner_log_mixt_density = log_omega;
    for (m in 1:M) {
    inner_log_mixt_density[m] += log_f[i, m];
    }
    log_weighted_mixt_density[k] += log_sum_exp(inner_log_mixt_density);
    log_weighted_mixt_density[k] += beta_proportion_lpdf(H[i] | mu[k], nu[k]);
    }
    target += log_sum_exp(log_weighted_mixt_density);
    }
    }
    "
    input_data <- list(n=sum(ind_est), K=K, M=M, comp_forecasts=component_forecasts[train_id,], comp_pits=component_pits[train_id,],
                       alpha_w=hyperpars$alpha_w,omega=hyperpars$omega_ew)
  }
  if((!grepl("EW",ensemble_name,fixed = TRUE)) && grepl("BMC1",ensemble_name,fixed = TRUE)){
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> K;
    int<lower=0> M;
    matrix[n,M] comp_forecasts;
    matrix[n,M] comp_pits;
  }
    
    transformed data {
    matrix[n,M] log_f = log(comp_forecasts);
    }
    
    parameters {
    vector<lower=0, upper=1>[K] mu;
    vector<lower=0>[K] nu;
    simplex[M] omega_k[K];
    simplex[K] w;
    }
    
    model {
    matrix[n, K] H;
    vector[K] log_weighted_mixt_density;
    vector[M] inner_log_mixt_density;
    vector[K] log_w = log(w);
    vector[M] log_omega_k[K];
    // priors
    for (k in 1:K) {
    log_omega_k[k] = log(omega_k[k]);
    H[, k] = comp_pits*omega_k[k];
    }
    for (i in 1:n) {
    log_weighted_mixt_density = log_w;
    for (k in 1:K) {
    inner_log_mixt_density = log_omega_k[k];
    for (m in 1:M) {
    inner_log_mixt_density[m] += log_f[i, m];
    }
    log_weighted_mixt_density[k] += log_sum_exp(inner_log_mixt_density);
    log_weighted_mixt_density[k] += beta_proportion_lpdf(H[i, k] | mu[k], nu[k]);
    }
    target += log_sum_exp(log_weighted_mixt_density);
    }
    }
    "
    input_data <- list(n=sum(ind_est), K=K, M=M, comp_forecasts=component_forecasts[train_id,], comp_pits=component_pits[train_id,])
  }
  if( grepl("EW",ensemble_name,fixed = TRUE) && grepl("BMC1",ensemble_name,fixed = TRUE)){
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> K;
    int<lower=0> M;
    vector<lower=0, upper=1>[M] omega;
    matrix[n,M] comp_forecasts;
    matrix[n,M] comp_pits;
  }
    
    transformed data {
    matrix[n,M] log_f = log(comp_forecasts);
    vector[M] log_omega = log(omega);
    vector[n] H = comp_pits*omega;
    }
    
    parameters {
    vector<lower=0, upper=1>[K] mu;
    vector<lower=0>[K] nu;
    simplex[K] w;
    }
    
    model {
    vector[K] log_weighted_mixt_density;
    vector[M] inner_log_mixt_density;
    vector[K] log_w = log(w);
    for (i in 1:n) {
    log_weighted_mixt_density = log_w;
    for (k in 1:K) {
    inner_log_mixt_density = log_omega;
    for (m in 1:M) {
    inner_log_mixt_density[m] += log_f[i, m];
    }
    log_weighted_mixt_density[k] += log_sum_exp(inner_log_mixt_density);
    log_weighted_mixt_density[k] += beta_proportion_lpdf(H[i] | mu[k], nu[k]);
    }
    target += log_sum_exp(log_weighted_mixt_density);
    }
    }
    "
    input_data <- list(n=sum(ind_est), K=K, M=M, comp_forecasts=component_forecasts[train_id,], comp_pits=component_pits[train_id,],
                       omega=hyperpars$omega_ew)
  }
  model_obj <- stan_model(model_code=model_type)
  pars <- rstan::optimizing(model_obj,data=input_data,seed=1234,hessian=T,init=initial_vals,
                            verbose=TRUE)
  
  if(grepl("EW-", ensemble_name, fixed = TRUE)){
    # ew-bmc
    omega <- rep(1/ncol(component_forecasts),ncol(component_forecasts))
    munu <- pars$par[1:(K*2)]
    w <- pars$par[(1+(2*K)):length(pars$par)]
    blp_train <- matrix(NA,ncol=K,nrow=nrow(component_forecasts[train_id,]))
    blp_test <- matrix(NA,ncol=K,nrow=nrow(component_forecasts[test_id,]))
    pbeta_mat_train <- matrix(NA,ncol=K,nrow=nrow(component_pits[train_id,]))
    pbeta_mat_test <- matrix(NA,ncol=K,nrow=nrow(component_pits[test_id,]))
    for (k in 1:K){
      assign(paste0("ab",k),c(munu[k]*munu[k+K],(1-munu[k])*munu[k+K]))
      h_train <- component_forecasts[train_id,] %*% omega
      h_test <- component_forecasts[test_id,] %*% omega
      beta_train <- dbeta(component_pits[train_id,] %*% omega,
                          shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
      beta_test <- dbeta(component_pits[test_id,] %*% omega,
                         shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
      blp_train[,k] <- h_train*beta_train
      blp_test[,k] <- h_test*beta_test
      pbeta_mat_train[,k] <- pbeta(component_pits[train_id,] %*% omega,
                             shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
      pbeta_mat_test[,k] <- pbeta(component_pits[test_id,] %*% omega,
                                   shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
    }
  } else {
    # BMC
    munu <- pars$par[1:(K*2)]
    omega <- pars$par[1+(K*2):(length(pars$par)-(K+1))]
    w <- pars$par[(1+(2*K)+(M*K)):length(pars$par)]
    blp_train <- matrix(NA,ncol=K,nrow=nrow(component_forecasts[train_id,]))
    blp_test <- matrix(NA,ncol=K,nrow=nrow(component_forecasts[test_id,]))
    pbeta_mat_train <- matrix(NA,ncol=K,nrow=nrow(component_pits[train_id,]))
    pbeta_mat_test <- matrix(NA,ncol=K,nrow=nrow(component_pits[test_id,]))
    for (k in 1:K){
      assign(paste0("ab",k),c(munu[k]*munu[k+K],(1-munu[k])*munu[k+K]))
      assign(paste0("omega",k),c())
      for (m in 1:M){
        assign(paste0("omega",k), c(get(paste0("omega",k)),omega[((m-1)*K)+k]))}
      h_train <- component_forecasts[train_id,] %*% c(get(paste0("omega",k)))
      h_test <- component_forecasts[test_id,] %*% c(get(paste0("omega",k)))
      beta_train <- dbeta(component_pits[train_id,] %*% get(paste0("omega",k)),
                    shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
      beta_test <- dbeta(component_pits[test_id,] %*% get(paste0("omega",k)),
                          shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
      blp_train[,k] <- h_train*beta_train
      blp_test[,k] <- h_test*beta_test
      pbeta_mat_train[,k] <- pbeta(component_pits[train_id,] %*% get(paste0("omega",k)),
                             shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
      pbeta_mat_test[,k] <- pbeta(component_pits[test_id,] %*% get(paste0("omega",k)),
                                   shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
    }
  }
  # make density (check google doc on idea of how to do it)
  density_train <- blp_train %*% w
  density_test <- blp_test %*% w
  ll_train     <- log(density_train)
  ll_test     <- log(density_test)
  PIT_train <- pbeta_mat_train %*% w
  PIT_test <- pbeta_mat_test %*% w
  results <- list(params=pars$par,train_ls=mean(ll_train),test_ls=mean(ll_test),train_pits=PIT_train,
                  test_pits=PIT_test, train_dens=density_train, test_dens = density_test)
  return(results)
}
  
##---------------------- Presentation functions (need to edit) --------------------------------##
make_ls <- function(forecast_names,scenario){
  mls <- matrix(NA, nrow = length(forecast_names), ncol = 2,
                dimnames = list(c(forecast_names),c("Training", "Test")))
  for(i in 1:length(forecast_names)){
    mls[i, 1] <- get(paste0(forecast_names[i],"_results_",scenario))$train_ls
    mls[i, 2] <- get(paste0(forecast_names[i],"_results_",scenario))$test_ls
    }
    mls<-round(mls, 3)
    return(mls)
}

make_table <- function(forecast_names,scenario,component_forecasts){
  table <- matrix(NA, nrow = 4, ncol = 5,
                  dimnames = list(c(forecast_names), c("$\\omega_1$", "$\\omega_2$", "$\\omega_3$", "$\\alpha$", "$\\beta$")))
  table[1,] <- c(get(paste0(forecast_names[1],"_results_",scenario))$params, rep(NA, 2))
  table[2,] <- c(get(paste0(forecast_names[2],"_results_",scenario))$params[-c(1:2)], 
                 get(paste0(forecast_names[2],"_results_",scenario))$params[1:2])
  table[3,] <- c(get(paste0(forecast_names[3],"_results_",scenario))$params, rep(NA, 2))
  table[4,] <- c(rep(1/ncol(component_forecasts),ncol(component_forecasts)), 
                 get(paste0(forecast_names[4],"_results_",scenario))$params)
  # table[5,] <- c(get(paste0(forecast_names[5],"_results_",scenario))$params[-c(1:2,6)], 
  #                get(paste0(forecast_names[5],"_results_",scenario))$params[1]*get(paste0(forecast_names[5],"_results_",scenario))$params[2],
  #                (1-get(paste0(forecast_names[5],"_results_",scenario))$params[1])*get(paste0(forecast_names[5],"_results_",scenario))$params[2])
  # table[6,] <- c(rep(1/ncol(component_forecasts),ncol(component_forecasts)), 
  #                get(paste0(forecast_names[6],"_results_",scenario))$params[1]*get(paste0(forecast_names[6],"_results_",scenario))$params[2],
  #                (1-get(paste0(forecast_names[6],"_results_",scenario))$params[1])*get(paste0(forecast_names[6],"_results_",scenario))$params[2])
  table<-round(table, 3)
  return(table)
}

make_table_bmc <- function(forecast_names,scenario){
  for(ens in forecast_names){
    if(grepl("EW-",ens,fixed=TRUE)){
      assign(paste0(gsub("EW-", "ew", ens),"_pars"),get(paste0(ens,"_results_",scenario))$params)
    } else {
      assign(paste0(ens,"_pars"),get(paste0(ens,"_results_",scenario))$params)
    }
  }
  
  table <- matrix(NA, nrow = 8, ncol = 30,
                  dimnames = list(c(forecast_names), c("$w_1$", "$w_2$", "$w_3$","$w_4$","$w_5$", "$\\alpha_1$", "$\\beta_1$",
                                                       "$\\alpha_2$", "$\\beta_2$", "$\\alpha_3$", "$\\beta_3$",
                                                       "$\\alpha_4$", "$\\beta_4$","$\\alpha_5$", "$\\beta_5$",
                                                       "$\\omega_{11}$", "$\\omega_{12}$", "$\\omega_{13}$",
                                                       "$\\omega_{21}$", "$\\omega_{22}$", "$\\omega_{23}$",
                                                       "$\\omega_{31}$", "$\\omega_{32}$", "$\\omega_{33}$",
                                                       "$\\omega_{41}$", "$\\omega_{42}$", "$\\omega_{43}$",
                                                       "$\\omega_{51}$", "$\\omega_{52}$", "$\\omega_{53}$")))
  table[1,] <- c(BMC2_pars[11:12],  rep(NA, 3), BMC2_pars[1]*BMC2_pars[3], (1-BMC2_pars[1])*BMC2_pars[3],
                 BMC2_pars[2]*BMC2_pars[4], (1-BMC2_pars[2])*BMC2_pars[4],rep(NA, 6),
                 BMC2_pars[5], BMC2_pars[7], BMC2_pars[9], BMC2_pars[6],BMC2_pars[8],BMC2_pars[10],rep(NA, 9))
  table[2,] <- c(ewBMC2_pars[5:6],  rep(NA, 3), ewBMC2_pars[1]*ewBMC2_pars[3], (1-ewBMC2_pars[1])*ewBMC2_pars[3],
                 ewBMC2_pars[2]*ewBMC2_pars[4], (1-ewBMC2_pars[2])*ewBMC2_pars[4],rep(NA, 6),rep(0.333, 6),rep(NA, 9))
  table[3,] <- c(BMC3_pars[16:18],  rep(NA, 2), BMC3_pars[1]*BMC3_pars[4], (1-BMC3_pars[1])*BMC3_pars[4],
                 BMC3_pars[2]*BMC3_pars[5], (1-BMC3_pars[2])*BMC3_pars[5],
                 BMC3_pars[3]*BMC3_pars[6], (1-BMC3_pars[3])*BMC3_pars[6],rep(NA, 4),
                 BMC3_pars[7], BMC3_pars[10], BMC3_pars[13], BMC3_pars[8],BMC3_pars[11],BMC3_pars[14],
                 BMC3_pars[9],BMC3_pars[12],BMC3_pars[15],rep(NA, 6))
  table[4,] <- c(ewBMC3_pars[7:9],  rep(NA, 2), ewBMC3_pars[1]*ewBMC3_pars[4], (1-ewBMC3_pars[1])*ewBMC3_pars[4],
                 ewBMC3_pars[2]*ewBMC3_pars[5], (1-ewBMC3_pars[2])*ewBMC3_pars[5],
                 ewBMC3_pars[3]*ewBMC3_pars[6], (1-ewBMC3_pars[3])*ewBMC3_pars[6],rep(NA, 4),rep(0.333, 9),rep(NA, 6))
  table[5,] <- c(BMC4_pars[21:24], rep(NA, 1), BMC4_pars[1]*BMC4_pars[5], (1-BMC4_pars[1])*BMC4_pars[5],
                 BMC4_pars[2]*BMC4_pars[6], (1-BMC4_pars[2])*BMC4_pars[6],
                 BMC4_pars[3]*BMC4_pars[7], (1-BMC4_pars[3])*BMC4_pars[7],
                 BMC4_pars[4]*BMC4_pars[8], (1-BMC4_pars[4])*BMC4_pars[8],rep(NA, 2),
                 BMC4_pars[9], BMC4_pars[13], BMC4_pars[17], BMC4_pars[10],BMC4_pars[14],BMC4_pars[18],
                 BMC4_pars[11],BMC4_pars[15],BMC4_pars[19],BMC4_pars[12],BMC4_pars[16],BMC4_pars[20],rep(NA, 3))
  table[6,] <- c(ewBMC4_pars[9:12],rep(NA, 1), ewBMC4_pars[1]*ewBMC4_pars[5], (1-ewBMC4_pars[1])*ewBMC4_pars[5],
                 ewBMC4_pars[2]*ewBMC4_pars[6], (1-ewBMC4_pars[2])*ewBMC4_pars[6],
                 ewBMC4_pars[3]*ewBMC4_pars[7], (1-ewBMC4_pars[3])*ewBMC4_pars[7],
                 ewBMC4_pars[4]*ewBMC4_pars[8], (1-ewBMC4_pars[4])*ewBMC4_pars[8],rep(NA, 2),
                 rep(0.333, 12),rep(NA, 3))
  table[7,] <- c(BMC5_pars[26:30], BMC5_pars[1]*BMC5_pars[6], (1-BMC5_pars[1])*BMC5_pars[6],
                 BMC5_pars[2]*BMC5_pars[7], (1-BMC5_pars[2])*BMC5_pars[7],
                 BMC5_pars[3]*BMC5_pars[8], (1-BMC5_pars[3])*BMC5_pars[8],
                 BMC5_pars[4]*BMC5_pars[9], (1-BMC5_pars[4])*BMC5_pars[9], 
                 BMC5_pars[5]*BMC5_pars[10], (1-BMC5_pars[5])*BMC5_pars[10], 
                 BMC5_pars[11], BMC5_pars[16], BMC5_pars[21], BMC5_pars[12],BMC5_pars[17],BMC5_pars[22],
                 BMC5_pars[13],BMC5_pars[18],BMC5_pars[23],BMC5_pars[14],BMC5_pars[19],BMC5_pars[24],
                 BMC5_pars[15],BMC5_pars[20],BMC5_pars[25])
  table[8,] <- c(ewBMC5_pars[11:15], ewBMC5_pars[1]*ewBMC5_pars[6], (1-ewBMC5_pars[1])*ewBMC5_pars[6],
                 ewBMC5_pars[2]*ewBMC5_pars[7], (1-ewBMC5_pars[2])*ewBMC5_pars[7],
                 ewBMC5_pars[3]*ewBMC5_pars[8], (1-ewBMC5_pars[3])*ewBMC5_pars[8],
                 ewBMC5_pars[4]*ewBMC5_pars[9], (1-ewBMC5_pars[4])*ewBMC5_pars[9],
                 ewBMC5_pars[5]*ewBMC5_pars[10], (1-ewBMC5_pars[5])*ewBMC5_pars[10],rep(0.333, 15))
  table<-round(table, 3)
  return(table)
} 

make_table2 <- function(forecast_names,scenario,component_forecasts){
  table <- matrix(NA, nrow = 4, ncol = 4,
                  dimnames = list(c(forecast_names), c("$\\omega_1$", "$\\omega_2$", "$\\alpha$", "$\\beta$")))
  table[1,] <- c(get(paste0(forecast_names[1],"_results_",scenario))$params, rep(NA, 2))
  table[2,] <- c(get(paste0(forecast_names[2],"_results_",scenario))$params[-c(1:2)], 
                 get(paste0(forecast_names[2],"_results_",scenario))$params[1:2])
  table[3,] <- c(get(paste0(forecast_names[3],"_results_",scenario))$params, rep(NA, 2))
  table[4,] <- c(rep(1/ncol(component_forecasts),ncol(component_forecasts)), 
                 get(paste0(forecast_names[4],"_results_",scenario))$params)
  # table[5,] <- c(get(paste0(forecast_names[5],"_results_",scenario))$params[-c(1:2,5)], 
  #                get(paste0(forecast_names[5],"_results_",scenario))$params[1]*get(paste0(forecast_names[5],"_results_",scenario))$params[2],
  #                (1-get(paste0(forecast_names[5],"_results_",scenario))$params[1])*get(paste0(forecast_names[5],"_results_",scenario))$params[2])
  # table[6,] <- c(rep(1/ncol(component_forecasts),ncol(component_forecasts)), 
  #                get(paste0(forecast_names[6],"_results_",scenario))$params[1]*get(paste0(forecast_names[6],"_results_",scenario))$params[2],
  #                (1-get(paste0(forecast_names[6],"_results_",scenario))$params[1])*get(paste0(forecast_names[6],"_results_",scenario))$params[2])
  table<-round(table, 3)
  return(table)
}

make_table_bmc2 <- function(forecast_names,scenario){
  for(ens in forecast_names){
    if(grepl("EW-",ens,fixed=TRUE)){
      assign(paste0(gsub("EW-", "ew", ens),"_pars"),get(paste0(ens,"_results_",scenario))$params)
    } else {
      assign(paste0(ens,"_pars"),get(paste0(ens,"_results_",scenario))$params)
    }
  }
  table <- matrix(NA, nrow = 8, ncol = 25,
                  dimnames = list(c(forecast_names), c("$w_1$", "$w_2$", "$w_3$","$w_4$","$w_5$", "$\\alpha_1$", "$\\beta_1$",
                                                       "$\\alpha_2$", "$\\beta_2$", "$\\alpha_3$", "$\\beta_3$",
                                                       "$\\alpha_4$", "$\\beta_4$","$\\alpha_5$", "$\\beta_5$",
                                                       "$\\omega_{11}$", "$\\omega_{12}$", 
                                                       "$\\omega_{21}$", "$\\omega_{22}$","$\\omega_{31}$", "$\\omega_{32}$",
                                                       "$\\omega_{41}$", "$\\omega_{42}$","$\\omega_{51}$", "$\\omega_{52}$")))
  table[1,] <- c(BMC2_pars[9:10],  rep(NA, 3), BMC2_pars[1]*BMC2_pars[3], (1-BMC2_pars[1])*BMC2_pars[3],
                 BMC2_pars[2]*BMC2_pars[4], (1-BMC2_pars[2])*BMC2_pars[4],rep(NA, 6),BMC2_pars[5], BMC2_pars[7],
                 BMC2_pars[6],BMC2_pars[8],rep(NA, 6))
  table[2,] <- c(ewBMC2_pars[5:6],  rep(NA, 3), ewBMC2_pars[1]*ewBMC2_pars[3], (1-ewBMC2_pars[1])*ewBMC2_pars[3],
                 ewBMC2_pars[2]*ewBMC2_pars[4], (1-ewBMC2_pars[2])*ewBMC2_pars[4],rep(NA, 6),rep(0.5, 4),rep(NA, 6))
  table[3,] <- c(BMC3_pars[13:15],  rep(NA, 2), BMC3_pars[1]*BMC3_pars[4], (1-BMC3_pars[1])*BMC3_pars[4],
                 BMC3_pars[2]*BMC3_pars[5], (1-BMC3_pars[2])*BMC3_pars[5],
                 BMC3_pars[3]*BMC3_pars[6], (1-BMC3_pars[3])*BMC3_pars[6],rep(NA, 4),
                 BMC3_pars[7], BMC3_pars[10], BMC3_pars[8],BMC3_pars[11], BMC3_pars[9],BMC3_pars[12],rep(NA, 4))
  table[4,] <- c(ewBMC3_pars[7:9],  rep(NA, 2), ewBMC3_pars[1]*ewBMC3_pars[4], (1-ewBMC3_pars[1])*ewBMC3_pars[4],
                 ewBMC3_pars[2]*ewBMC3_pars[5], (1-ewBMC3_pars[2])*ewBMC3_pars[5],
                 ewBMC3_pars[3]*ewBMC3_pars[6], (1-ewBMC3_pars[3])*ewBMC3_pars[6],rep(NA, 4),rep(0.5, 6),rep(NA, 4))
  table[5,] <- c(BMC4_pars[17:20],rep(NA, 1), BMC4_pars[1]*BMC4_pars[5], (1-BMC4_pars[1])*BMC4_pars[5],
                 BMC4_pars[2]*BMC4_pars[6], (1-BMC4_pars[2])*BMC4_pars[6],
                 BMC4_pars[3]*BMC4_pars[7], (1-BMC4_pars[3])*BMC4_pars[7],
                 BMC4_pars[4]*BMC4_pars[8], (1-BMC4_pars[4])*BMC4_pars[8],rep(NA, 2),
                 BMC4_pars[9], BMC4_pars[13], BMC4_pars[10],BMC4_pars[14],BMC4_pars[11],
                 BMC4_pars[15],BMC4_pars[12],BMC4_pars[16],rep(NA, 2))
  table[6,] <- c(ewBMC4_pars[9:12],rep(NA, 1), 
                 ewBMC4_pars[1]*ewBMC4_pars[5], (1-ewBMC4_pars[1])*ewBMC4_pars[5],
                 ewBMC4_pars[2]*ewBMC4_pars[6], (1-ewBMC4_pars[2])*ewBMC4_pars[6],
                 ewBMC4_pars[3]*ewBMC4_pars[7], (1-ewBMC4_pars[3])*ewBMC4_pars[7],
                 ewBMC4_pars[4]*ewBMC4_pars[8], (1-ewBMC4_pars[4])*ewBMC4_pars[8],
                 rep(NA, 2),rep(0.5, 8),rep(NA, 2))
  table[7,] <- c(BMC4_pars[21:25], BMC4_pars[1]*BMC4_pars[6], (1-BMC4_pars[1])*BMC4_pars[6],
                 BMC4_pars[2]*BMC4_pars[7], (1-BMC4_pars[2])*BMC4_pars[7],
                 BMC4_pars[3]*BMC4_pars[8], (1-BMC4_pars[3])*BMC4_pars[8],
                 BMC4_pars[4]*BMC4_pars[9], (1-BMC4_pars[4])*BMC4_pars[9],
                 BMC4_pars[5]*BMC4_pars[10], (1-BMC4_pars[5])*BMC4_pars[10],
                 BMC4_pars[11], BMC4_pars[16], BMC4_pars[12],BMC4_pars[17],BMC4_pars[13],
                 BMC4_pars[18],BMC4_pars[14],BMC4_pars[19],BMC4_pars[15],BMC4_pars[20])
  table[8,] <- c(ewBMC4_pars[11:15], 
                 ewBMC4_pars[1]*ewBMC4_pars[6], (1-ewBMC4_pars[1])*ewBMC4_pars[6],
                 ewBMC4_pars[2]*ewBMC4_pars[7], (1-ewBMC4_pars[2])*ewBMC4_pars[7],
                 ewBMC4_pars[3]*ewBMC4_pars[8], (1-ewBMC4_pars[3])*ewBMC4_pars[8],
                 ewBMC4_pars[4]*ewBMC4_pars[9], (1-ewBMC4_pars[4])*ewBMC4_pars[9],
                 ewBMC4_pars[5]*ewBMC4_pars[10], (1-ewBMC4_pars[5])*ewBMC4_pars[10],
                 rep(0.5, 10))

  table<-round(table, 3)
  return(table)
}

plot_PITs <- function(comp_probs,forecast_names,scenario,train=TRUE){
  par(mfrow = c(5, 3),mai = c(0.1, 0.5, 0.5, 0.1),mai = c(0.2, 0.3, 0.3, 0.2))
  if(train){
  # component model 
  for (i in 1:ncol(comp_probs)){
    hist(comp_probs[ind_est,i], breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 2), main = paste0("train_component ",i),xlab="",cex.main=0.8)
    abline( h = 1, lty = 2)
  }
  # ensembles
  for (ens in forecast_names){
    hist(get(paste0(ens,"_results_",scenario))$train_pits, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 2), main = paste0("train_",ens),xlab="",cex.main=0.8)
    abline( h = 1, lty = 2)
    
  }} else{
    for (i in 1:ncol(comp_probs)){
      hist(comp_probs[!ind_est,i], breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
           ylim = c(0, 2), main = paste0("test_component ",i),xlab="",cex.main=0.8)
      abline( h = 1, lty = 2)
    }
    # ensembles
    for (ens in forecast_names){
      hist(get(paste0(ens,"_results_",scenario))$test_pits, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
           ylim = c(0, 2), main = paste0("test_",ens),xlab="",cex.main=0.8)
      abline( h = 1, lty = 2)
    }
  }
}

plot_PDFs <- function(Y,comp_dens,forecast_names,scenario, yaxis){
  par(mfrow = c(5, 3),mai = c(0.1, 0.5, 0.5, 0.1),mai = c(0.2, 0.3, 0.3, 0.2))
  # reorder Y based on low to high here --- doesn't quite work
  ord_y <- Y[!ind_est][order(Y[!ind_est])]
  # get index
      for (i in 1:ncol(comp_dens)){
        plot(ord_y,comp_dens[!ind_est,i][order(Y[!ind_est])], type="l",
             main = paste0("test_component",i),xlab="",cex.main=0.8, ylim=yaxis)
      }
      # ensembles
      for (ens in forecast_names){
        plot(ord_y,get(paste0(ens,"_results_",scenario))$test_dens[order(Y[!ind_est])], type="l",
             main = paste0("test_",ens),xlab="",cex.main=0.8, ylim=yaxis)
    
    }
}

plot_PDFs_s1 <- function(Y,comp_dens,forecast_names,scenario){
  par(mfrow = c(5, 3),mai = c(0.1, 0.5, 0.5, 0.1),mai = c(0.2, 0.3, 0.3, 0.2))
  # reorder Y based on low to high here --- doesn't quite work
  ord_y <- Y[!ind_est][order(Y[!ind_est])]
  # get index
  for (i in 1:ncol(comp_dens)){
    plot(ord_y,comp_dens[!ind_est,i][order(Y[!ind_est])], type="l",
         main = paste0("test_component",i),xlab="",cex.main=0.8, col='grey')
    lines(predict(loess(comp_dens[!ind_est,i][order(Y[!ind_est])]~ord_y)), lwd=2)
  }
  # ensembles
  for (ens in forecast_names){
    plot(ord_y,get(paste0(ens,"_results_",scenario))$test_dens[order(Y[!ind_est])], type="l",
         main = paste0("test_",ens),xlab="",cex.main=0.8, col='grey')
    lines(predict(loess(get(paste0(ens,"_results_",scenario))$test_dens[order(Y[!ind_est])]~ord_y)), lwd=2)
  }
}