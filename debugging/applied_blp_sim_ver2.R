#---------------------------------------------- input functions ------------------------------------------------#
cdf_to_pdf <- function(cdf_vector){
  prob_vector <- c(cdf_vector[1],c(cdf_vector - dplyr::lag(cdf_vector))[-1])
  pdf <- prob_vector/sum(prob_vector)
  return(pdf)
}

#-------------------------------------------------- ensemble functions -----------------------------------------------#
# get cdf and params from here and make pdf 
# after make pdf, use pdf to calculate ls and pit outside for both train and test
make_ensemble <- function(ensemble_name, K, M, train_indices, truths, component_pdfs, 
                          initial_vals, hyperpars=NULL){
  # generate some inputs
  # make cdfs
  component_cdfs <- component_pdfs %>%
    dplyr::group_by(model,index) %>%
    dplyr::mutate(cdf_val=cumsum(binned_probs)) %>%
    ungroup()
  # get likelihood for training 
  train_lik <- component_pdfs %>%
    dplyr::left_join(Y, by="index") %>%
    dplyr::filter(index %in% train_indices) %>%
    dplyr::group_by(model,index) %>%
    dplyr::filter((Y>=bin_start) & (Y< bin_end)) %>%
    ungroup() %>%
    dplyr::select("model","binned_probs","index") %>%
    tidyr::pivot_wider(id_cols=index,
                       names_from=model,
                       values_from = binned_probs) %>%
    dplyr::select(-"index") %>%
    as.matrix()
  # get pit for training
  dat <- component_cdfs %>%
    dplyr::left_join(Y, by="index") %>%
    dplyr::filter(index %in% train_indices)
  train_pit <- dat %>%
    dplyr::group_by(model,index) %>%
    dplyr::filter((Y>=bin_start) & (Y< bin_end)) %>%
    ungroup() %>%
    dplyr::select("model","cdf_val","index") %>%
    tidyr::pivot_wider(id_cols=index,
                       names_from=model,
                       values_from = cdf_val) %>%
    dplyr::select(-"index") %>%
    as.matrix()
  train_prepit <- dat %>%
    dplyr::group_by(model,index) %>%
    dplyr::filter((Y>=bin_start) & (Y > bin_end)) %>%
    dplyr::filter(bin_start==max(bin_start)) %>% 
    ungroup() %>%
    dplyr::select("model","cdf_val","index") %>%
    tidyr::pivot_wider(id_cols=index,
                       names_from=model,
                       values_from = cdf_val) %>%
    dplyr::select(-"index") %>%
    as.matrix()
  if(ensemble_name=="TLP") {
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> M;
    matrix[n,M] comp_likelihood;
    }
    parameters {
    simplex[M] omega;
    }
    
    model {
    vector[M] log_omega = log(omega);
    matrix[n, M] log_likelihood = log(comp_likelihood);
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
    input_data <- list(n=length(train_indices), M=M, comp_likelihood=train_lik)
  }
  else if(ensemble_name=="BLP"){
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> K;
    int<lower=0> M;
    matrix[n,M] comp_pits;
    matrix[n,M] comp_prepits;
    }

    parameters {
    vector<lower=0, upper=1>[K] mu;
    vector<lower=0>[K] nu;
    simplex[M] omega_k[K];
    simplex[K] w;
    }
 
    model {
    matrix[n, K] H_a;
    matrix[n, K] H_b;
    vector[K] log_weighted_mixt_cdf_a;
    vector[K] log_weighted_mixt_cdf_b;
    vector[K] log_w = log(w);
    
    for (k in 1:K) {
    H_a[, k] = comp_pits*omega_k[k];
    H_b[, k] = comp_prepits*omega_k[k];
    }
    for (i in 1:n) {
    log_weighted_mixt_cdf_a = log_w;
    log_weighted_mixt_cdf_b = log_w;
    for (k in 1:K) {
    log_weighted_mixt_cdf_a[k] += beta_proportion_lcdf(H_a[i,k] | mu[k], nu[k]);
    log_weighted_mixt_cdf_b[k] += beta_proportion_lcdf(H_b[i,k] | mu[k], nu[k]);
    }
    target += log_diff_exp(log_sum_exp(log_weighted_mixt_cdf_a), log_sum_exp(log_weighted_mixt_cdf_b));
    }
    }
    "
    input_data <- list(n=length(train_indices), K=1, M=M, comp_pits=train_pit, comp_prepits=train_prepit)
  }
  else if(ensemble_name=="EW-BLP"){
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> K;
    int<lower=0> M;
    vector<lower=0, upper=1>[M] omega;
    matrix[n,M] comp_pits;
    matrix[n,M] comp_prepits;
    }
    
    transformed data {
    vector[n] H_a = comp_pits*omega;
    vector[n] H_b = comp_prepits*omega;
    }
    
    parameters {
    vector<lower=0, upper=1>[K] mu;
    vector<lower=0>[K] nu;
    simplex[K] w;
    }
    
    model {
    vector[K] log_weighted_mixt_cdf_a;
    vector[K] log_weighted_mixt_cdf_b;
    vector[K] log_w = log(w);
    
    for (i in 1:n) {
    log_weighted_mixt_cdf_a = log_w;
    log_weighted_mixt_cdf_b = log_w;
    for (k in 1:K) {
    log_weighted_mixt_cdf_a[k] += beta_proportion_lcdf(H_a[i] | mu[k], nu[k]);
    log_weighted_mixt_cdf_b[k] += beta_proportion_lcdf(H_b[i] | mu[k], nu[k]);
    }
    target += log_diff_exp(log_sum_exp(log_weighted_mixt_cdf_a), log_sum_exp(log_weighted_mixt_cdf_b));
    }
    }
    "
    input_data <- list(n=length(train_indices), K=1, M=M, omega=rep(1/M,M),
                       comp_pits=train_pit, comp_prepits=train_prepit)
    
  }
  else if(grepl("BMC",ensemble_name,fixed = TRUE) && !grepl("EW-BMC",ensemble_name,fixed = TRUE)){
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> K;
    int<lower=0> M;
    vector<lower=0>[K] alpha_w;
    vector<lower=0>[M] alpha_omega;
    matrix[n,M] comp_pits;
    matrix[n,M] comp_prepits;
    }

    parameters {
    vector<lower=0, upper=1>[K] mu;
    vector<lower=0>[K] nu;
    simplex[M] omega_k[K];
    simplex[K] w;
    }

    model {
    matrix[n, K] H_a;
    matrix[n, K] H_b;
    vector[K] log_weighted_mixt_cdf_a;
    vector[K] log_weighted_mixt_cdf_b;
    vector[K] log_w = log(w);
    // priors
    for (k in 1:K) {
    mu[k] ~ beta(2, 2);
    nu[k] ~ gamma(0.1, 0.1);
    omega_k[k] ~ dirichlet(alpha_omega);
    }
    w ~ dirichlet(alpha_w);

    for (k in 1:K) {
    H_a[, k] = comp_pits*omega_k[k];
    H_b[, k] = comp_prepits*omega_k[k];
    }
    for (i in 1:n) {
    log_weighted_mixt_cdf_a = log_w;
    log_weighted_mixt_cdf_b = log_w;
    for (k in 1:K) {
    log_weighted_mixt_cdf_a[k] += beta_proportion_lcdf(H_a[i, k] | mu[k], nu[k]);
    log_weighted_mixt_cdf_b[k] += beta_proportion_lcdf(H_b[i, k] | mu[k], nu[k]);
    }
    target += log_diff_exp(log_sum_exp(log_weighted_mixt_cdf_a), log_sum_exp(log_weighted_mixt_cdf_b));
    }
    }
    "
    input_data <- list(n=length(train_indices), K=K, M=M,
                       comp_pits=train_pit,comp_prepits=train_prepit,
                       alpha_w=hyperpars$alpha_w,alpha_omega=hyperpars$alpha_omega)

    }
  else if(grepl("EW-BMC",ensemble_name,fixed = TRUE)) {
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> K;
    int<lower=0> M;
    vector<lower=0, upper=1>[M] omega;
    matrix[n,M] comp_pits;
    matrix[n,M] comp_prepits;
    }
    transformed data {
    vector[n] H_a = comp_pits*omega;
    vector[n] H_b = comp_prepits*omega;
    }
    parameters {
    vector<lower=0, upper=1>[K] mu;
    vector<lower=0>[K] nu;
    simplex[K] w;
    }
      
    model {
    vector[K] log_w = log(w);
    // priors
    for (k in 1:K) {
    mu[k] ~ beta(2, 2);
    nu[k] ~ gamma(0.1, 0.1);
    }
    w ~ dirichlet(alpha_w);
    for (i in 1:n) {
    log_weighted_mixt_cdf_a = log_w;
    log_weighted_mixt_cdf_b = log_w;
    for (k in 1:K) {
    log_weighted_mixt_cdf_a[k] += beta_proportion_lcdf(H_a[i] | mu[k], nu[k]);
    log_weighted_mixt_cdf_b[k] += beta_proportion_lcdf(H_b[i] | mu[k], nu[k]);
    }
    target += log_diff_exp(log_sum_exp(log_weighted_mixt_cdf_a), log_sum_exp(log_weighted_mixt_cdf_b));
    }
    }
    "
    input_data <- list(n=length(train_indices), K=K, M=M,
                       comp_pits=train_pit,
                       comp_prepits=train_prepit,
                       alpha_w=hyperpars$alpha_w,omega=hyperpars$omega_ew)
  }
  model_obj <- stan_model(model_code=model_type)
  pars <- rstan::optimizing(model_obj,data=input_data,seed=1234,hessian=T,init=initial_vals,
                            verbose=TRUE)
  #make a list of results
  params <- pars$par
  return(params)
}

# build ensemble pdfs
make_pdfs <- function(ensemble_name, params, K,component_pdfs){
  # generate some inputs
  M <- unique(component_pdfs$model)
  component_cdfs <- component_pdfs %>%
    dplyr::group_by(model,index) %>%
    dplyr::mutate(cdf_val=cumsum(binned_probs)) %>%
    ungroup() %>%
    dplyr::select(-"binned_probs")
  comp_matrix <- component_pdfs %>%
    tidyr::pivot_wider(id_cols =c(bin_start,bin_end,index), 
                       names_from = model, 
                       values_from = binned_probs) %>%
    as.matrix()
  comp_cdf_matrix <- component_cdfs %>%
    dplyr::select("bin_start","bin_end","model","cdf_val","index") %>%
    tidyr::pivot_wider(id_cols =c(bin_start,bin_end,index),
                       names_from = model,
                       values_from = binned_probs) %>%
    as.matrix()
  # frame for pdf
  ensemble_pdf <- data.frame(comp_matrix[,1:4])
  # build train and test cdf/pdf part
  if(ensemble_name=="EW"|ensemble_name=="TLP"){
    pdf_vals <- comp_matrix[,4:(4+M)] %*% params
    ensemble_pdf$binned_probs <- pdf_vals
    ensemble_pdf <- ensemble_pdf %>%
      dplyr::group_by(index) %>%
      dplyr::mutate(binned_probs=binned_probs/sum(binned_probs)) %>%
      ungroup()
  } 
  else if(grepl("EW-", ensemble_name, fixed = TRUE)){
    # ew-bmc or ew-blp
    ## extract params
    omega <- rep(1/M,M)
    munu <- params[1:(K*2)]
    w <- params[(1+(2*K)):length(params)]
    # equally-weighed component models (ew ensemble)
    H <- comp_cdf_matrix[,4:(4+M)] %*% omega
    # build matrices for each k
    pbeta_mat <- matrix(NA,ncol=K,nrow=nrow(comp_cdf_matrix))
    for (k in 1:K){
      assign(paste0("ab",k),c(munu[k]*munu[k+K],(1-munu[k])*munu[k+K]))
      pbeta_mat[,k] <- pbeta( H, shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
    }
    # turn pbeta_mat into a combined mixture (pbeta_mat %*% w)
    ensemble_pdf$cdf_vals <- pbeta_mat %*% w
    # then turn  a combined mixture into pdf using the cdf to pdf function
    ensemble_pdf <- ensemble_pdf %>%
      dplyr::group_by(index) %>%
      dplyr::mutate(binned_probs=cdf_to_pdf(cdf_vals)) %>%
      ungroup() %>%
      dplyr::select(-"cdf_vals")
  } else {
    # BMC or blp
    ## extract params
    munu <- params[1:(K*2)]
    omega <- params[1+(K*2):(length(params)-(K+1))]
    w <- params[(1+(2*K)+(M*K)):length(params)]
    # build matrices for each k
    pbeta_mat <- matrix(NA,ncol=K,nrow=nrow(comp_cdf_matrix))
    for (k in 1:K){
      assign(paste0("ab",k),c(munu[k]*munu[k+K],(1-munu[k])*munu[k+K]))
      assign(paste0("omega",k),c())
      for (m in 1:M){
        assign(paste0("omega",k), c(get(paste0("omega",k)),omega[((m-1)*K)+k]))
      }
      H <- comp_cdf_matrix[,4:(4+M)] %*% c(get(paste0("omega",k)))
      pbeta_mat[,k] <- pbeta( H, shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
      }
    # turn pbeta_mat into a combined mixture (pbeta_mat %*% w)
    ensemble_pdf$cdf_vals <- pbeta_mat %*% w
    # then turn  a combined mixture into pdf using the cdf to pdf function
    ensemble_pdf <- ensemble_pdf %>%
      dplyr::group_by(index) %>%
      dplyr::mutate(binned_probs=cdf_to_pdf(cdf_vals)) %>%
      ungroup() %>%
      dplyr::select(-"cdf_vals")
  }
  return(ensemble_pdf)
}

# a function to calculate ls and pits
ls_pit_calc <- function(pdf, Y, train_indices){
  pit_frame <- pdf %>%
    dplyr::group_by(index) %>%
    dplyr::mutate(pitv=cumsum(binned_probs)) %>%
    ungroup() %>%
    dplyr::left_join(Y, by="index") %>%
    dplyr::filter((Y>=bin_start) & (Y< bin_end)) 
  ls_frame <- pdf %>%
    dplyr::left_join(Y, by="index") %>%
    dplyr::group_by(index) %>%
    dplyr::filter((Y>=bin_start) & (Y< bin_end)) %>%
    ungroup() %>%
    dplyr::mutate(ls=ifelse(binned_probs==0,-10,log(binned_probs)))
  # join
  frame <- pit_frame %>%
    dplyr::left_join(ls_frame, by=c("bin_start","bin_end","index")) 
  score_frame_train <- frame  %>%
    dplyr::filter(index %in% train_indices) %>%
    dplyr::mutate(mean_ls=mean(ls))
  score_frame_test <- frame  %>%
    dplyr::filter(!(index %in% train_indices)) %>%
    dplyr::mutate(mean_ls=mean(ls))
  results <- list(mean_ls_train=unique(score_frame_train$mean_ls),
                  mean_ls_test=unique(score_frame_test$mean_ls),
                  pit_train=score_frame_train$pitv,
                  pit_test=score_frame_test$pitv)
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


#### debugging area
## try re-writing this to identity function (normal beta dis)
make_ensemble2 <- function(ensemble_name, K, M, train_indices, truths, component_pdfs,
                          initial_vals, hyperpars=NULL){
  # generate some inputs
  # make cdfs
  component_cdfs <- component_pdfs %>%
    dplyr::group_by(model,index) %>%
    dplyr::mutate(cdf_val=cumsum(binned_probs)) %>%
    ungroup()
  # get likelihood for training
  train_lik <- component_pdfs %>%
    dplyr::left_join(Y, by="index") %>%
    dplyr::filter(index %in% train_indices) %>%
    dplyr::group_by(model,index) %>%
    dplyr::filter((Y>=bin_start) & (Y< bin_end)) %>%
    ungroup() %>%
    dplyr::select("model","binned_probs","index") %>%
    tidyr::pivot_wider(id_cols=index,
                       names_from=model,
                       values_from = binned_probs) %>%
    dplyr::select(-"index") %>%
    as.matrix()
  # get pit for training
  dat <- component_cdfs %>%
    dplyr::left_join(Y, by="index") %>%
    dplyr::filter(index %in% train_indices)
  train_pit <- dat %>%
    dplyr::group_by(model,index) %>%
    dplyr::filter((Y>=bin_start) & (Y< bin_end)) %>%
    ungroup() %>%
    dplyr::select("model","cdf_val","index") %>%
    tidyr::pivot_wider(id_cols=index,
                       names_from=model,
                       values_from = cdf_val) %>%
    dplyr::select(-"index") %>%
    as.matrix()
  train_prepit <- dat %>%
    dplyr::group_by(model,index) %>%
    dplyr::filter((Y>=bin_start) & (Y > bin_end)) %>%
    dplyr::filter(bin_start==max(bin_start)) %>%
    ungroup() %>%
    dplyr::select("model","cdf_val","index") %>%
    tidyr::pivot_wider(id_cols=index,
                       names_from=model,
                       values_from = cdf_val) %>%
    dplyr::select(-"index") %>%
    as.matrix()
  if(ensemble_name=="BLP"){
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> K;
    int<lower=0> M;
    matrix[n,M] comp_pits;
    matrix[n,M] comp_prepits;
    }

    parameters {
    vector<lower=0>[K] alpha;
    vector<lower=0>[K] beta;
    simplex[M] omega_k[K];
    simplex[K] w;
    }
    
    model {
    matrix[n, K] H_a;
    matrix[n, K] H_b;
    vector[K] log_weighted_mixt_cdf_a;
    vector[K] log_weighted_mixt_cdf_b;
    vector[K] log_w = log(w);

    for (k in 1:K) {
    H_a[, k] = comp_pits*omega_k[k];
    H_b[, k] = comp_prepits*omega_k[k];
    }
    for (i in 1:n) {
    log_weighted_mixt_cdf_a = log_w;
    log_weighted_mixt_cdf_b = log_w;
    for (k in 1:K) {
    log_weighted_mixt_cdf_a[k] += beta_lcdf(H_a[i, k] | alpha[k], beta[k]);
    log_weighted_mixt_cdf_b[k] += beta_lcdf(H_b[i, k] | alpha[k], beta[k]);
    }
    target += log_diff_exp(log_sum_exp(log_weighted_mixt_cdf_a), log_sum_exp(log_weighted_mixt_cdf_b));
    }
    }
    "
    input_data <- list(n=length(train_indices), M=M, comp_pits=train_pit, comp_prepits=train_prepit)
  }
  else if(ensemble_name=="EW-BLP"){
    model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> M;
    vector<lower=0, upper=1>[M] omega;
    matrix[n,M] comp_pits;
    matrix[n,M] comp_prepits;
    }
    transformed data {
    vector[n] H_a = comp_pits*omega;
    vector[n] H_b = comp_prepits*omega;
    }
    parameters {
    real<lower=0> alpha;
    real<lower=0> beta;
    }

    model {
    for (i in 1:n) {
    target += log_diff_exp(beta_lcdf(H_a[i] | alpha, beta),beta_lcdf(H_b[i] | alpha, beta));
    }
    }
    "
    input_data <- list(n=length(train_indices), M=M, omega=rep(1/M,M),
                       comp_pits=train_pit, comp_prepits=train_prepit)

  }
  
  model_obj <- stan_model(model_code=model_type,verbose = T,)
  pars <- rstan::optimizing(model_obj,data=input_data,seed=1234, hessian = T,init=initial_vals,
                            verbose=TRUE)
  # make a list of results
  return(params)
}