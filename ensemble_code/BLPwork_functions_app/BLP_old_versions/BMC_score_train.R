library(tidyverse)
library(rstan)
library(foreach)
library(doParallel)

# get mixture params on the whole training seasons
# info
args <- commandArgs(TRUE)
test_s <- as.numeric(args[1])
if(test_s == 17){
  train_season<-c("2010/2011","2011/2012","2012/2013","2013/2014","2014/2015","2015/2016")
} else if (test_s == 18){
  train_season<-c("2010/2011","2011/2012","2012/2013","2013/2014","2014/2015","2015/2016","2016/2017")
} else if (test_s == 19){
  train_season<-c("2010/2011","2011/2012","2012/2013","2013/2014","2014/2015","2015/2016","2016/2017","2017/2018")
}
targets<-c("1 wk ahead","2 wk ahead","3 wk ahead","4 wk ahead")
# set path
path <- "~/ch1dat/"
#path <- "/Users/ahctun_woon/git/end-of-2019-2020/cdc-flusight-ensemble/"
load(file = paste0("BLP",test_s,".rda"))
load(file = paste0("EW_BLP",test_s,".rda"))
# get in scores and pits
pitByModel<-read.csv(paste0(path,"pit_modelcol.csv")) %>%
  dplyr::filter(calendar_week %in% c(43:53,1:18),
                season %in% train_season)
prepitByModel<-read.csv(paste0(path,"prepit_modelcol.csv")) %>%
  dplyr::filter(calendar_week %in% c(43:53,1:18),
                season %in% train_season)
# pitByModel<-read.csv(paste0(path,"calibration_work/BLPdata_app/pit_modelcol.csv")) %>%
#   dplyr::filter(calendar_week %in% c(43:53,1:18),
#                 season %in% train_season)
# prepitByModel<-read.csv(paste0(path,"calibration_work/BLPdata_app/prepit_modelcol.csv")) %>%
#   dplyr::filter(calendar_week %in% c(43:53,1:18),
#                 season %in% train_season)
# info
componentModel_list <- colnames(pitByModel)[5:ncol(pitByModel)]
mod_num <- length(componentModel_list)

##--------------------------------------------------Functions to get params-------------------------------------------##
make_ensemble_emp <- function(ensemble_name, K, M, component_pits, component_prepits){
  # generate some inputs
  # get pit for training
  train_pit <- component_pits[,5:ncol(component_pits)] %>%
    as.matrix()
  train_prepit <- component_prepits[,5:ncol(component_prepits)] %>%
    as.matrix()
  n_obs <- nrow(train_pit)
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
    input_data <- list(n=n_obs, K=K, M=M,comp_pits=train_pit,comp_prepits=train_prepit)
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
    input_data <- list(n=n_obs, K=K, M=M, omega=rep(1/M,M),comp_pits=train_pit,comp_prepits=train_prepit)
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

# cv function 
cv_func <- function(tar,validation_season,ensemble_name,K,M,component_pits,component_prepits){
  train_cv <- setdiff(train_season,validation_season)
  component_pits_cv <- component_pits %>%
    dplyr::filter(season %in% train_cv,
                  target == tar)
  component_prepits_cv <- component_prepits %>%
    dplyr::filter(season %in% train_cv,
                  target == tar)
  pars <- make_ensemble_emp(ensemble_name, K=K, M=M, 
                            component_pits=component_pits_cv, 
                            component_prepits=component_prepits_cv)
  return(pars)
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
  if(grepl("EW_", ensemble_name, fixed = TRUE)){
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
##----------------------------------------------------- BLP -----------------------------------------------------------##

#k=1
EW_BLP_params <- lapply(1:4, function(x) cv_func(targets[x],
                                                 NULL,
                                                 "EW-BLP",
                                                 K=1,
                                                 M=mod_num,
                                                 component_pits=pitByModel,
                                                 component_prepits=prepitByModel))


# cv 
for (i in 2:5){
  for (j in 1:4){
    assign(paste0("EW_BMC",i,"_params_cv_target",j),
           lapply(train_season, 
                  function(x) cv_func(targets[j],
                                      x,
                                      "EW-BLP",
                                      K=i,
                                      M=mod_num,
                                      component_pits=pitByModel,
                                      component_prepits=prepitByModel)
           )
    )
  }
}
##---------------------------------------------- save outputs ---------------------------------------------------------##
assign(paste0("EW_BLP",test_s),list(EW_BLP_params=EW_BLP_params,
                                    EW_BMC2_params_cv_target1=EW_BMC2_params_cv_target1,
                                    EW_BMC3_params_cv_target1=EW_BMC3_params_cv_target1,
                                    EW_BMC4_params_cv_target1=EW_BMC4_params_cv_target1,
                                    EW_BMC5_params_cv_target1=EW_BMC5_params_cv_target1,
                                    EW_BMC2_params_cv_target2=EW_BMC2_params_cv_target2,
                                    EW_BMC3_params_cv_target2=EW_BMC3_params_cv_target2,
                                    EW_BMC4_params_cv_target2=EW_BMC4_params_cv_target2,
                                    EW_BMC5_params_cv_target2=EW_BMC5_params_cv_target2,
                                    EW_BMC2_params_cv_target3=EW_BMC2_params_cv_target3,
                                    EW_BMC3_params_cv_target3=EW_BMC3_params_cv_target3,
                                    EW_BMC4_params_cv_target3=EW_BMC4_params_cv_target3,
                                    EW_BMC5_params_cv_target3=EW_BMC5_params_cv_target3,
                                    EW_BMC2_params_cv_target4=EW_BMC2_params_cv_target4,
                                    EW_BMC3_params_cv_target4=EW_BMC3_params_cv_target4,
                                    EW_BMC4_params_cv_target4=EW_BMC4_params_cv_target4,
                                    EW_BMC5_params_cv_target4=EW_BMC5_params_cv_target4))
save(get(paste0("EW_BLP",test_s)),file=paste0("EW_BLP",test_s,".rda"))