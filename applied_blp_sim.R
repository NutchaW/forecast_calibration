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

# make mixture ensembles
make_ensemble <- function(ensemble_name,K=NULL,component_forecasts,component_pits, initial_vals=NULL,
                          train_id, test_id){
  M <- ncol(component_forecasts)
  n_obs <- nrow(component_forecasts[train_id,])
  if(ensemble_name=="TLP") {
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> M;
    matrix[n,M] comp_forecasts;
    }
    parameters {
    simplex[M] omega_k;
    }
    
    model {
    vector[M] log_omega = log(omega_k);
    matrix[n, M] log_likelihood = log(comp_forecasts);
    vector[M] log_weighted_mixt_pdf;
    for (i in 1:n) {
    log_weighted_mixt_pdf = log_omega;
    for (m in 1:M) {
    log_weighted_mixt_pdf[m] += log_likelihood[i, m];
    }
    target += log_sum_exp(log_weighted_mixt_pdf);
    }
    }
    "
    input_data <- list(n=n_obs, M=M, comp_forecasts=component_forecasts[train_id,])
  }
  else if((!grepl("EW_",ensemble_name,fixed = TRUE)) && grepl("BMC",ensemble_name,fixed = TRUE)){
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
    vector<lower=0>[K] alpha;
    vector<lower=0>[K] beta;
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
    log_weighted_mixt_density[k] += beta_lpdf(H[i, k] | alpha[k], beta[k]);
    }
    target += log_sum_exp(log_weighted_mixt_density);
    }
    }
    "
    input_data <- list(n=n_obs, K=K, M=M, comp_forecasts=component_forecasts[train_id,], comp_pits=component_pits[train_id,])
  }
  else if(grepl("EW_",ensemble_name,fixed = TRUE)&& grepl("BMC",ensemble_name,fixed = TRUE)){
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
    vector<lower=0>[K] alpha;
    vector<lower=0>[K] beta;
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
    log_weighted_mixt_density[k] += beta_lpdf(H[i] | alpha[k], beta[k]);
    }
    target += log_sum_exp(log_weighted_mixt_density);
    }
    }
    "
    input_data <- list(n=n_obs, K=K, M=M, comp_forecasts=component_forecasts[train_id,], 
                       comp_pits=component_pits[train_id,],
                       omega=rep(1/M,M))
  }
  if(ensemble_name!="EW"){
  model_obj <- stan_model(model_code=model_type)
  pars <- rstan::optimizing(model_obj,data=input_data,seed=1234,hessian=T,init=initial_vals,verbose=TRUE)
  }
  if(ensemble_name=="TLP"|ensemble_name=="EW"){
    if(ensemble_name=="TLP"){
      omega <- pars$par
    } 
    else {
      omega <- rep(1/M,M)
    }
    density_train <- component_forecasts[train_id,] %*% omega
    density_test <- component_forecasts[test_id,] %*% omega
    PIT_train <- component_pits[train_id,] %*% omega
    PIT_test <- component_pits[test_id,] %*% omega
  }
  else {
    if(!grepl("EW_", ensemble_name, fixed = TRUE)){
      omega <- pars$par[1+(K*2):(length(pars$par)-(K+1))]
      w <- pars$par[(1+(2*K)+(M*K)):length(pars$par)]
    }
    else{
      omega <- rep(1/M,M*K)
      w <- pars$par[(1+(2*K)):length(pars$par)]
    }
    ab <- pars$par[1:(K*2)]
    blp_train <- matrix(NA,ncol=K,nrow=nrow(component_forecasts[train_id,]))
    blp_test <- matrix(NA,ncol=K,nrow=nrow(component_forecasts[test_id,]))
    pbeta_mat_train <- matrix(NA,ncol=K,nrow=nrow(component_pits[train_id,]))
    pbeta_mat_test <- matrix(NA,ncol=K,nrow=nrow(component_pits[test_id,]))
    for (k in 1:K){
      ab_k <- c(ab[k],ab[k+K])
      omega_k <- c()
      for (m in 1:M){
        omega_k <- c(omega_k,omega[((m-1)*K)+k])
        }
      h_train <- component_forecasts[train_id,] %*% omega_k
      h_test <- component_forecasts[test_id,] %*% omega_k
      beta_train <- dbeta(component_pits[train_id,] %*% omega_k,shape1 = ab_k[1], shape2 = ab_k[2])
      beta_test <- dbeta(component_pits[test_id,] %*% omega_k,shape1 = ab_k[1], shape2 = ab_k[2])
      blp_train[,k] <- h_train*beta_train
      blp_test[,k] <- h_test*beta_test
      pbeta_mat_train[,k] <- pbeta(component_pits[train_id,] %*% omega_k,
                                   shape1 = ab_k[1], shape2 = ab_k[2])
      pbeta_mat_test[,k] <- pbeta(component_pits[test_id,] %*% omega_k,
                                  shape1 = ab_k[1], shape2 = ab_k[2])
    }
    density_train <- blp_train %*% w
    density_test <- blp_test %*% w
    PIT_train <- pbeta_mat_train %*% w
    PIT_test <- pbeta_mat_test %*% w
  }
  ll_train     <- log(density_train)
  ll_test     <- log(density_test)
  if(ensemble_name=="EW"){
  results <- list(params=omega,train_ls=mean(ll_train),test_ls=mean(ll_test),train_pits=PIT_train,
                  test_pits=PIT_test, train_dens=density_train, test_dens = density_test)
  } 
  else {
  results <- list(params=pars$par,train_ls=mean(ll_train),test_ls=mean(ll_test),train_pits=PIT_train,
                  test_pits=PIT_test, train_dens=density_train, test_dens = density_test)
  }
  return(results)
}

# make a wrapper for running 
non_mixture_wrapper <- function(method, component_forecasts, component_pits,train_indices){
  omega_vec_inits <- sample.int(100, size = ncol(component_forecasts), replace = FALSE)
  omega <- omega_vec_inits/sum(omega_vec_inits)
  if(method=="TLP"|method=="EW"){
    if(method=="TLP"){
      init_vals <- list(omega_k=omega)
      }
    else{
      init_vals <- NULL
      }
    k_param <- NULL
    }
  else {
    init_k1 <- list(alpha=rnorm(1,1,1e-2),beta=rnorm(1,1,1e-2),w=1.0)
    dim(init_k1$alpha) <- 1
    dim(init_k1$beta) <- 1
    dim(init_k1$w) <- 1
    k_param <- 1
    if(method=="EW_BMC1"){
      init_vals <- init_k1
      } 
    else {
      init_vals <- c(init_k1,list(omega_k=omega))
      dim(init_vals$omega_k) <- c(1,ncol(component_forecasts))
    }
    }
  res_list <- make_ensemble(method,K=k_param,component_forecasts, component_pits, 
                            init_vals,train_indices, !train_indices)
  return(res_list)
}

cv_wrapper <- function(method_list,cv_list,component_forecasts,component_pits){
  omega_vec_inits <- sample.int(100, size = ncol(component_forecasts), replace = FALSE)
  omega <- omega_vec_inits/sum(omega_vec_inits)
  # do 5-cv
  cv_info <- data.frame(cbind(method=rep(method_list,length(cv_list)),
                              fold=sort(rep(1:length(cv_list),length(method_list))), 
                              train_ls=NA, validation_ls=NA))
  for(k in 1:length(cv_list)){
    for(mix_ens in method_list){
      K <- as.numeric(substr(mix_ens, nchar(mix_ens), nchar(mix_ens)))
      w_vec_inits <- sample.int(100, size = K, replace = FALSE)
      init_vals <- list(alpha=rnorm(K,1,1e-2),beta=rnorm(K,1,1e-2),w=w_vec_inits/sum(w_vec_inits))
      if(!grepl("EW_", mix_ens, fixed = TRUE)){
        init_vals <- c(init_vals,list(omega_k=t(matrix(rep(omega,K),ncol=K))))
        }
    score <- make_ensemble(mix_ens,K,component_forecasts,component_pits, init_vals,
                               c(1:(n-(n/5)))[-cv_list[[k]]],cv_list[[k]])
    cv_info$train_ls[cv_info$method== mix_ens & cv_info$fold== k] <- score$train_ls
    cv_info$validation_ls[cv_info$method== mix_ens & cv_info$fold== k] <- score$test_ls
    }
    }
  cv_info <- cv_info %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(mean_valid_ls=mean(as.numeric(validation_ls)),
                  mean_train_ls=mean(as.numeric(train_ls)))  %>%
    dplyr::ungroup() 
  # split the data set into ew and bmc and find ones with max valid ls
  # max valid within 1 sd
  ewbmc <- cv_info[grepl("EW_", cv_info$method, fixed = TRUE),] %>%
    dplyr::mutate(lower_bound=max(mean_valid_ls)-sd(mean_valid_ls),
                  num_k = as.numeric(substr(method, nchar(method), nchar(method))),
                  within_sd = ifelse(mean_valid_ls>=lower_bound,TRUE,FALSE)) %>%
    dplyr::select(-c("train_ls","validation_ls","fold")) %>%
    distinct() %>%
    dplyr::filter(within_sd==TRUE)
  bmc <- cv_info[!grepl("EW_", cv_info$method, fixed = TRUE),] %>%
    dplyr::mutate(lower_bound=max(mean_valid_ls)-sd(mean_valid_ls),
                  num_k = as.numeric(substr(method, nchar(method), nchar(method))),
                  within_sd = ifelse(mean_valid_ls>=lower_bound,TRUE,FALSE)) %>%
    dplyr::select(-c("train_ls","validation_ls","fold")) %>%
    distinct() %>%
    dplyr::filter(within_sd==TRUE)
  method_ewbmc <- ewbmc[which.min(ewbmc$num_k),]$method
  method_bmc <- bmc[which.min(bmc$num_k),]$method
  return(list(select=c(method_bmc,method_ewbmc),cv_info=cv_info))
}


mixture_wrapper <- function(method, component_forecasts, component_pits,train_indices){
  omega_vec_inits <- sample.int(100, size = ncol(component_forecasts), replace = FALSE)
  omega <- omega_vec_inits/sum(omega_vec_inits)
  K <- as.numeric(substr(method, nchar(method), nchar(method)))
  w_vec_inits <- sample.int(100, size = K, replace = FALSE)
  init_vals <- list(alpha=rnorm(K,1,1e-2),beta=rnorm(K,1,1e-2),w=w_vec_inits/sum(w_vec_inits))
  if(!grepl("EW_", method, fixed = TRUE)){
    init_vals <- c(init_vals,list(omega_k=t(matrix(rep(omega,K),ncol=K))))
    } 
  res_list <- make_ensemble(method,K,component_forecasts,component_pits, init_vals,
                            train_indices,!train_indices)
  return(res_list)
}
##---------------------- Presentation functions (need to edit) --------------------------------##
make_ls <- function(forecast_names,scenario){
  name_list <- replace(forecast_names, forecast_names %in% c("EW_BMC1","BMC1"), c("BLP","EW_BLP"))
  name_list <- gsub("_", "-", name_list)
  mls <- matrix(NA, nrow = length(forecast_names), ncol = 2,
                dimnames = list(c(name_list),c("Training", "Test")))
  for(i in 1:length(forecast_names)){
    mls[i, 1] <- get(paste0(forecast_names[i],"_results_",scenario))$train_ls
    mls[i, 2] <- get(paste0(forecast_names[i],"_results_",scenario))$test_ls
    }
  mls<-round(mls, 3)
  return(mls)
}

make_table <- function(forecast_names,scenario){
  # make table and kable it
  param_names <- c("$w_1$", "$w_2$", "$w_3$","$w_4$","$w_5$", 
      "$\\alpha_1$", "$\\beta_1$", "$\\alpha_2$", "$\\beta_2$", "$\\alpha_3$", "$\\beta_3$", 
      "$\\alpha_4$", "$\\beta_4$","$\\alpha_5$", "$\\beta_5$",
      "$\\omega_{11}$", "$\\omega_{12}$", "$\\omega_{13}$",
      "$\\omega_{21}$", "$\\omega_{22}$", "$\\omega_{23}$",
      "$\\omega_{31}$", "$\\omega_{32}$", "$\\omega_{33}$",
      "$\\omega_{41}$", "$\\omega_{42}$", "$\\omega_{43}$",
      "$\\omega_{51}$", "$\\omega_{52}$", "$\\omega_{53}$")
  param_cols <- c("w[1]", "w[2]", "w[3]","w[4]","w[5]", 
                  "alpha[1]", "beta[1]", "alpha[2]", "beta[2]", "alpha[3]", "beta[3]", 
                  "alpha[4]", "beta[4]","alpha[5]", "beta[5]",
                  "omega_k[1,1]", "omega_k[1,2]", "omega_k[1,3]",
                  "omega_k[2,1]", "omega_k[2,2]", "omega_k[2,3]",
                  "omega_k[3,1]", "omega_k[3,2]", "omega_k[3,3]",
                  "omega_k[4,1]", "omega_k[4,2]", "omega_k[4,3]",
                  "omega_k[5,1]", "omega_k[5,2]", "omega_k[5,3]")
  table_params <- data.frame(matrix(NA, nrow = 6, ncol = length(param_cols)))
  colnames(table_params) <- param_cols
  table_params$Method <- forecast_names
  table_params <- table_params[c(ncol(table_params),1:(ncol(table_params)-1))]
  for(ens in forecast_names){
    params <- round(get(paste0(ens,"_results_s",scenario))$params,3)
    for(i in 1:length(params)){
      if(ens=="TLP"|ens=="EW"){
        table_params[table_params$Method==ens, paste0("omega_k[",1,",",i,"]")] <- params[i]
      }
      else{
        table_params[table_params$Method==ens, names(params)[i]] <- params[i]
      }
      if(grepl("EW_", ens, fixed = TRUE)){
        K <- as.numeric(substr(ens, nchar(ens), nchar(ens)))
        M<- ifelse(scenario!=3,3,2)
        for(k in 1:K){
          for(m in 1:M){
          table_params[table_params$Method==ens, paste0("omega_k[",k,",",m,"]")] <- round(1/M,3)
          }
        }
      }
      }
  }
  table_params[,2][table_params$Method=="EW_BMC1"|table_params$Method=="BMC1"] <- NA
  table_params$Method[table_params$Method=="EW_BMC1"] <- "EW_BLP"
  table_params$Method[table_params$Method=="BMC1"] <- "BLP"
  names(table_params)[2:ncol(table_params)] <- param_names
  table_params <- table_params[rowSums(is.na(table_params))<ncol(table_params),
                               colSums(is.na(table_params))<nrow(table_params)]
  table_params$Method <- gsub("_", "-", table_params$Method)
  return(table_params)
} 


plot_PITs <- function(comp_probs,forecast_names,scenario,train=TRUE){
  name_list <- replace(forecast_names, forecast_names %in% c("EW_BMC1","BMC1"), c("BLP","EW_BLP"))
  par(mfrow = c(3, 3))
  if(scenario!="s1"){
    if(train){
      # component model 
      for (i in 1:ncol(comp_probs)){
        hist(comp_probs[ind_est,i], breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
             ylim = c(0, 4), main = paste0("train_component ",i),xlab="",cex.main=0.8)
        abline( h = 1, lty = 2)
      }
      # ensembles
      for (j in 1:length(forecast_names)){
        hist(get(paste0(forecast_names[j],"_results_",scenario))$train_pits, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
             ylim = c(0, 4), main = paste0(name_list[j]),xlab="",cex.main=0.8)
        abline( h = 1, lty = 2)
        
      }} else{
        for (i in 1:ncol(comp_probs)){
          hist(comp_probs[!ind_est,i], breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
               ylim = c(0, 4), main = paste0("component ",i),xlab="",cex.main=0.8)
          abline( h = 1, lty = 2)
        }
        # ensembles
        for (j in 1:length(forecast_names)){
          hist(get(paste0(forecast_names[j],"_results_",scenario))$test_pits, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
               ylim = c(0, 4), main = paste0(name_list[j]),xlab="",cex.main=0.8)
          abline( h = 1, lty = 2)
        }
      }
  }
  else{
    if(train){
      # component model 
      for (i in 1:ncol(comp_probs)){
        hist(comp_probs[ind_est,i], breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
             ylim = c(0, 2), main = paste0("train_component ",i),xlab="",cex.main=0.8)
        abline( h = 1, lty = 2)
      }
      # ensembles
      for (j in 1:length(forecast_names)){
        hist(get(paste0(forecast_names[j],"_results_",scenario))$train_pits, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
             ylim = c(0, 2), main = paste0("train_",name_list[j]),xlab="",cex.main=0.8)
        abline( h = 1, lty = 2)
        
      }} else{
        for (i in 1:ncol(comp_probs)){
          hist(comp_probs[!ind_est,i], breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
               ylim = c(0, 2), main = paste0("component ",i),xlab="",cex.main=0.8)
          abline( h = 1, lty = 2)
        }
        # ensembles
        for (j in 1:length(forecast_names)){
          hist(get(paste0(forecast_names[j],"_results_",scenario))$test_pits, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
               ylim = c(0, 2), main = paste0(name_list[j]),xlab="",cex.main=0.8)
          abline( h = 1, lty = 2)
        }
      }
  }
}

plot_PDFs <- function(Y,comp_dens,forecast_names,scenario, yaxis){
   name_list <- replace(forecast_names, forecast_names %in% c("EW_BMC1","BMC1"), c("BLP","EW_BLP"))
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
   name_list <- replace(forecast_names, forecast_names %in% c("EW_BMC1","BMC1"), c("BLP","EW_BLP"))
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

# maybe?? make plot beta cdf function

make_table2 <- function(forecast_names, res_list){
  # make table and kable it
  param_names <- c("$w_1$", "$w_2$", "$w_3$","$w_4$","$w_5$", 
                   "$\\alpha_1$", "$\\beta_1$", "$\\alpha_2$", "$\\beta_2$", "$\\alpha_3$", "$\\beta_3$", 
                   "$\\alpha_4$", "$\\beta_4$","$\\alpha_5$", "$\\beta_5$",
                   "$\\omega_{11}$", "$\\omega_{12}$", "$\\omega_{13}$",
                   "$\\omega_{21}$", "$\\omega_{22}$", "$\\omega_{23}$",
                   "$\\omega_{31}$", "$\\omega_{32}$", "$\\omega_{33}$",
                   "$\\omega_{41}$", "$\\omega_{42}$", "$\\omega_{43}$",
                   "$\\omega_{51}$", "$\\omega_{52}$", "$\\omega_{53}$")
  param_cols <- c("w[1]", "w[2]", "w[3]","w[4]","w[5]", 
                  "alpha[1]", "beta[1]", "alpha[2]", "beta[2]", "alpha[3]", "beta[3]", 
                  "alpha[4]", "beta[4]","alpha[5]", "beta[5]",
                  "omega_k[1,1]", "omega_k[1,2]", "omega_k[1,3]",
                  "omega_k[2,1]", "omega_k[2,2]", "omega_k[2,3]",
                  "omega_k[3,1]", "omega_k[3,2]", "omega_k[3,3]",
                  "omega_k[4,1]", "omega_k[4,2]", "omega_k[4,3]",
                  "omega_k[5,1]", "omega_k[5,2]", "omega_k[5,3]")
  table_params <- data.frame(matrix(NA, nrow = 6, ncol = length(param_cols)))
  colnames(table_params) <- param_cols
  table_params$Method <- forecast_names
  table_params <- table_params[c(ncol(table_params),1:(ncol(table_params)-1))]
  for(ens in forecast_names){
    params <- round(res_list[[ens]],3)
    for(i in 1:length(params)){
      if(ens=="TLP"|ens=="EW"){
        table_params[table_params$Method==ens, paste0("omega_k[",1,",",i,"]")] <- params[i]
      }
      else{
        table_params[table_params$Method==ens, names(params)[i]] <- params[i]
      }
      if(grepl("EW_", ens, fixed = TRUE)){
        K <- as.numeric(substr(ens, nchar(ens), nchar(ens)))
        M<- 3
        for(k in 1:K){
          for(m in 1:M){
            table_params[table_params$Method==ens, paste0("omega_k[",k,",",m,"]")] <- round(1/M,3)
          }
        }
      }
    }
  }
  table_params[,2][table_params$Method=="EW_BMC1"|table_params$Method=="BMC1"] <- NA
  table_params$Method[table_params$Method=="EW_BMC1"] <- "EW_BLP"
  table_params$Method[table_params$Method=="BMC1"] <- "BLP"
  names(table_params)[2:ncol(table_params)] <- param_names
  table_params <- table_params[rowSums(is.na(table_params))<ncol(table_params),
                               colSums(is.na(table_params))<nrow(table_params)]
  table_params$Method <- gsub("_", "-", table_params$Method)
  return(table_params)
} 


#### plot ensemble inside and beta mixture before combine

plot_details <- function(Y,table_row){
  ord_y <- Y[!ind_est][order(Y[!ind_est])]
  # blp
  omegas <- sapply(1:3, function(x) table_row[2,13+x])
  inside_ens <- comp_forecasts2[!ind_est,][order(Y[!ind_est]),] %*% omegas
  ## ew-blp
  ab <- sapply(1:2, function(x) table_row[4,5+x])
  inside_ens2 <- comp_forecasts2[!ind_est,][order(Y[!ind_est]),] %*% rep(1/3,3)
 # plot
  par(mfrow = c(2, 2))
  # blp
  # plot(ord_y,inside_ens[,1], type="l",
  #      main = "Ensemble inside beta transformation",xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
  plot(ord_y,get(paste0("BMC1_results_s2"))$test_dens[order(Y[!ind_est])], type="l",
       main = "BLP-Ensemble inside beta transformation",xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
  plot(ord_y,get(paste0("BMC1_results_s2"))$test_dens[order(Y[!ind_est])], type="l",
       main = paste0("BLP-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
  plot(ord_y,inside_ens2[,1], type="l",
       main = "EW_BLP-Ensemble inside beta transformation",xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
  plot(ord_y,get(paste0("EW_BMC1_results_s2"))$test_dens[order(Y[!ind_est])], type="l",
       main = paste0("EW_BLP-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
}

# plot_details2 <- function(Y,table_row){
#   ord_y <- Y[!ind_est][order(Y[!ind_est])]
#   pdf <- comp_forecasts2[!ind_est,][order(Y[!ind_est]),]
#   ## bmc3
#   ab1 <- sapply(1:2, function(x) table_row[5,5+x])
#   ab2 <- sapply(3:4, function(x) table_row[5,5+x])
#   # ab3 <- sapply(5:6, function(x) table_row[5,5+x])
#   omegas1 <- sapply(1:3, function(x) table_row[5,13+x])
#   omegas2 <- sapply(4:6, function(x) table_row[5,13+x])
#   # omegas3 <- sapply(5:6, function(x) table_row[5,13+x])
#   inside_ens1 <- as.vector(pdf %*% omegas1)
#   inside_ens2 <- as.vector(pdf %*% omegas2)
#   # inside_ens3 <- pdf %*% omegas3
#   beta_ens1 <- dbeta(pdf %*% omegas1,shape1=ab1[1],shape2=ab1[2])
#   beta_ens2 <- dbeta(pdf %*% omegas2,shape1=ab2[1],shape2=ab2[2])
#   # beta_ens3 <- dbeta(pdf %*% omegas3,shape1=ab3[1],shape2=ab3[2])
#   # ew-bmc4
#   ab11 <- sapply(1:2, function(x) table_row[6,5+x])
#   ab21 <- sapply(3:4, function(x) table_row[6,5+x])
#   ab31 <- sapply(5:6, function(x) table_row[6,5+x])
#   ab41 <- sapply(7:8, function(x) table_row[6,5+x])
#   ab51 <- sapply(9:10, function(x) table_row[6,5+x])
#   inside_ens11 <- pdf %*% rep(1/3,3)
#   beta_ens11 <- dbeta(inside_ens11,shape1=ab11[1],shape2=ab11[2])
#   beta_ens21 <- dbeta(inside_ens11,shape1=ab21[1],shape2=ab21[2])
#   beta_ens31 <- dbeta(inside_ens11,shape1=ab31[1],shape2=ab31[2])
#   beta_ens41 <- dbeta(inside_ens11,shape1=ab41[1],shape2=ab41[2])
#   beta_ens51 <- dbeta(inside_ens11,shape1=ab51[1],shape2=ab51[2])
#   # plot
#   par(mfrow = c(5, 2),byrow=T)
#   # bmc3
#   plot(ord_y,inside_ens1, type="l",
#        main = "BMC2-Ensemble inside beta transformation",xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
#   plot(ord_y,inside_ens2, type="l",
#        main = "BMC2-Ensemble inside beta transformation",xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
#   # plot(ord_y,inside_ens3[,1], type="l",
#   #      main = "BMC3-Ensemble inside beta transformation",xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
#   plot(ord_y,beta_ens1, type="l",
#        main = paste0("BMC2-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
#   plot(ord_y,beta_ens2, type="l",
#        main = paste0("BMC2-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
#   # plot(ord_y,beta_ens3, type="l",
#   #      main = paste0("BMC3-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,1.5))
#   #ew-bmc4
#   plot(ord_y,inside_ens11[,1], type="l",
#        main = "EW_BMC5-Ensemble inside beta transformation",xlab="",ylab="",cex.main=0.8, ylim=c(0,1))
#   plot(ord_y,beta_ens11 , type="l",
#        main = paste0("EW_BMC5-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,5))
#   plot(ord_y,beta_ens21 , type="l",
#        main = paste0("EW_BMC5-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,1))
#   plot(ord_y,beta_ens31 , type="l",
#        main = paste0("EW_BMC5-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,1))
#   plot(ord_y,beta_ens41 , type="l",
#        main = paste0("EW_BMC5-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,5))
#   plot(ord_y,beta_ens51 , type="l",
#        main = paste0("EW_BMC5-Beta-transformed ensemble"),xlab="",ylab="",cex.main=0.8, ylim=c(0,5))
# }
