#---------------------------------------------- input functions ------------------------------------------------#
cdf_to_pdf <- function(cdf_vector){
  prob_vector <- c(cdf_vector[1],c(cdf_vector - dplyr::lag(cdf_vector))[-1])
  pdf <- prob_vector/sum(prob_vector)
  return(pdf)
}

#-------------------------------------------------- ensemble functions -----------------------------------------------#
make_ensemble_emp <- function(ensemble_name, K, M, Y, train_indices, component_pdfs, 
                              hyperpars=NULL, test_indices){
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
    w_vec_inits <- sample.int(100, size = 3, replace = FALSE)
    init_vals <- list(omega=c(w_vec_inits/sum(w_vec_inits)))
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
      vector<lower=0.0001,upper=0.999>[K] alpha;
      vector<lower=0.0001,upper=65>[K] beta;
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
          print(\"i = \", i);
          print(\"Hs: \");
          print(H_a[i, k]);
          print(H_b[i, k]);
          print(\"Hs_truncated: \");
          print((H_a[i, k]*0.9943) + 0.00015);
          print((H_b[i, k]*0.9943) + 0.00015);
          print(alpha[k]);
          print(beta[k]);
          log_weighted_mixt_cdf_a[k] += beta_proportion_lcdf((H_a[i, k] * 0.9943) + 0.00015| alpha[k], beta[k]);
          log_weighted_mixt_cdf_b[k] += beta_proportion_lcdf((H_b[i, k] * 0.9943) + 0.00015| alpha[k], beta[k]);
        }
        target += log_diff_exp(log_sum_exp(log_weighted_mixt_cdf_a), log_sum_exp(log_weighted_mixt_cdf_b));
      }
    }
    "
    w_vec_inits <- sample.int(100, size = K, replace = FALSE)
    init_vals <- list(
      # alpha=rnorm(K,1,1e-3),
      # beta=rnorm(K,1,1e-3),
      alpha=rnorm(K,0.5,1e-3),
      beta=rnorm(K,2,1e-3),
      w=w_vec_inits/sum(w_vec_inits)
      )
    if(K==1){
      om_inits <- sample.int(10, size = M, replace = FALSE)
      dim(init_vals$alpha) <- 1
      dim(init_vals$beta) <- 1
      dim(init_vals$w) <- 1
      init_vals$omega_k <- om_inits/sum(om_inits)
      dim(init_vals$omega_k) <- c(1,M)
    } 
    else{
      om_inits <- sample.int(10, size = M, replace = FALSE)
      init_vals$omega_k <- t(replicate(K, om_inits/sum(om_inits)))
    }
    input_data <- list(n=length(train_indices), K=K, M=M,
                       comp_pits=train_pit,comp_prepits=train_prepit)
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
      vector<lower=0.0001,upper=0.999>[K] alpha;
      vector<lower=0.0001,upper=65>[K] beta;
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
          print(\"i = \", i);
          print(\"Hs: \");
          print(H_a[i]);
          print(H_b[i]);
          print(\"Hs_truncated: \");
          print((H_a[i]*0.9943) + 0.00015);
          print((H_b[i]*0.9943) + 0.00015);
          print(alpha[k]);
          print(beta[k]);
          log_weighted_mixt_cdf_a[k] += beta_proportion_lcdf((H_a[i] * 0.9943) + 0.00015| alpha[k], beta[k]);
          log_weighted_mixt_cdf_b[k] += beta_proportion_lcdf((H_b[i] * 0.9943) + 0.00015| alpha[k], beta[k]);
        }
        target += log_diff_exp(log_sum_exp(log_weighted_mixt_cdf_a), log_sum_exp(log_weighted_mixt_cdf_b));
      }
    }
    "
    w_vec_inits <- sample.int(100, size = K, replace = FALSE)
    init_vals <- list(
      # alpha=rnorm(K,1,1e-3),
      # beta=rnorm(K,1,1e-3),
      alpha=rnorm(K,0.5,1e-3),
      beta=rnorm(K,2,1e-3),
      w=w_vec_inits/sum(w_vec_inits),
      alpha=rnorm(K,0.5,1e-3),
      beta=rnorm(K,2,1e-3)
      )
    if(K==1){
    dim(init_vals$alpha) <- K
    dim(init_vals$beta) <- K
    dim(init_vals$w) <- K
    }
    input_data <- list(n=length(train_indices), K=K, M=M, omega=rep(1/M,M),
                       comp_pits=train_pit,comp_prepits=train_prepit)
    
  }

  model_obj <- stan_model(model_code=model_type)
  pars <- rstan::optimizing(model_obj,data=input_data,hessian=T,
                            init=init_vals, 
                            seed=1234,
                            verbose=TRUE)
  #make a list of results
  params <- pars$par
  return(params)
}

# build ensemble pdfs
make_pdfs <- function(ensemble_name, params, K,component_pdfs,M){
  if(ensemble_name=="BLP"|ensemble_name=="EW_BLP"){
    sub("BMC1","BLP",ensemble_name)
  }
  # generate some inputs
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
                       values_from = cdf_val) %>%
    as.matrix()
  # frame for pdf
  ensemble_pdf <- data.frame(comp_matrix[,1:3])
  # build train and test cdf/pdf part
  if(ensemble_name=="EW"|ensemble_name=="TLP"){
    # 3 is number of information columns
    pdf_vals <- comp_matrix[,4:ncol(comp_matrix)] %*% params
    ensemble_pdf$binned_probs <- pdf_vals
    ensemble_pdf <- ensemble_pdf %>%
      dplyr::group_by(index) %>%
      dplyr::mutate(binned_probs=binned_probs/sum(binned_probs)) %>%
      ungroup()
  }
  else if(grepl("EW_", ensemble_name, fixed = TRUE)){
    # ew-bmc or ew-blp
    ## extract params
    omega <- rep(1/M,M)
    munu <- params[1:(K*2)]
    w <- params[(1+(2*K)):length(params)]
    # equally-weighed component models (ew ensemble)
    H <- comp_cdf_matrix[,4:ncol(comp_matrix)] %*% omega
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
      H <- comp_cdf_matrix[,4:ncol(comp_cdf_matrix)] %*% c(get(paste0("omega",k)))
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

# build cv frame

cv_frame <- function(names,component_pdfs,train_indices,valid_index_list,M,Y){
  cv_binned <- data.frame(cbind(method=names,
                                 mean_train_ls=rep(NA,length(names)),
                                 mean_valid_ls=rep(NA,length(names))
  )
  )
  for(m_ens in c("BM","EW_BM")){
    if(m_ens=="BM"){
      k <- c(1,3,5,7)
    }else{
      k <- c(2,4,6,8)
    }
    for(i in 2:5){
      # for each fold
      ls_train <- c()
      ls_test <- c()
      for(f in 1:5){
        pdf <- make_pdfs(m_ens, get(paste0(m_ens,i,"_params"))[[f]], i,component_pdfs,M=M)
        res <- ls_pit_calc(pdf, Y, train_indices[-c(valid_index_list[[f]])])
        ls_train <- c(ls_train,res$mean_ls_train)
        ls_test <- c(ls_test,res$mean_ls_test)
      }
      cv_binned[k[i-1],2] <- mean(ls_train)
      cv_binned[k[i-1],3] <- mean(ls_test)
    }
  } 
  cv_binned$mean_train_ls <- as.numeric(cv_binned$mean_train_ls )
  cv_binned$mean_valid_ls <- as.numeric(cv_binned$mean_valid_ls)
  ewbmc <- cv_binned[grepl("EW_", cv_binned$method, fixed = TRUE),] %>%
    dplyr::mutate(lower_bound=(max(mean_valid_ls)-sd(mean_valid_ls)),
                  num_k = as.numeric(substr(method, nchar(method), nchar(method))),
                  within_sd = ifelse(mean_valid_ls>=lower_bound,TRUE,FALSE)) %>%
    dplyr::filter(within_sd==TRUE)
  bmc <- cv_binned[!grepl("EW_", cv_binned$method, fixed = TRUE),] %>%
    dplyr::mutate(lower_bound=(max(mean_valid_ls)-sd(mean_valid_ls)),
                  num_k = as.numeric(substr(method, nchar(method), nchar(method))),
                  within_sd = ifelse(mean_valid_ls>=lower_bound,TRUE,FALSE)) %>%
    dplyr::filter(within_sd==TRUE)
  method_ewbmc <- ewbmc[which.min(ewbmc$num_k),]$method
  method_bmc <- bmc[which.min(bmc$num_k),]$method
  return(list(select=c(method_bmc,method_ewbmc),cv_info=cv_binned))
}

make_ls_bin <- function(forecast_names){
  name_list1 <- replace(forecast_names, forecast_names %in% c("EW_BMC1","BMC1"), c("BLP","EW_BLP"))
  name_list <- gsub("_", "-", name_list1)
  mls <- matrix(NA, nrow = length(forecast_names), ncol = 2,
                dimnames = list(c(name_list),c("Training", "Test")))
  for(i in 1:length(name_list1)){
    mls[i, 1] <- get(paste0(name_list1[i],"_results_binned"))$mean_ls_train
    mls[i, 2] <- get(paste0(name_list1[i],"_results_binned"))$mean_ls_test
  }
  mls<-round(mls, 3)
  return(mls)
}

make_table_bin <- function(forecast_names,M){
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
    if(ens=="BMC1"){
      params <- round(get(paste0("BLP_params")),3)
    } else if (ens=="EW_BMC1"){
      params <- round(get(paste0("EW_BLP_params")),3)
    } else{
      params <- round(get(paste0(ens,"_params")),3)
    }
    for(i in 1:length(params)){
      if(ens=="TLP"|ens=="EW"){
        table_params[table_params$Method==ens, paste0("omega_k[",1,",",i,"]")] <- params[i]
      }
      else{
        table_params[table_params$Method==ens, names(params)[i]] <- params[i]
      }
      if(grepl("EW_", ens, fixed = TRUE)){
        K <- as.numeric(substr(ens, nchar(ens), nchar(ens)))
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

