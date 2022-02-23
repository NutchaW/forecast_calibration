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
          //print(\"i = \", i);
          //print(\"Hs: \");
          //print(H_a[i, k]);
          //print(H_b[i, k]);
          //print(\"Hs_truncated: \");
          //print((H_a[i, k]*0.9943) + 0.00015);
          //print((H_b[i, k]*0.9943) + 0.00015);
          //print(alpha[k]);
          //print(beta[k]);
          log_weighted_mixt_cdf_a[k] += beta_proportion_lcdf((H_a[i, k] * 0.9943) + 0.00015| alpha[k], beta[k]);
          log_weighted_mixt_cdf_b[k] += beta_proportion_lcdf((H_b[i, k] * 0.9943) + 0.00015| alpha[k], beta[k]);
        }
        target += log_diff_exp(log_sum_exp(log_weighted_mixt_cdf_a), log_sum_exp(log_weighted_mixt_cdf_b));
      }
    }
    "
    w_vec_inits <- sample.int(100, size = K, replace = TRUE)
    init_vals <- list(
      alpha=rnorm(K,0.5,1e-3),
      beta=rnorm(K,2,1e-3),
      w=w_vec_inits/sum(w_vec_inits)
    )
    if(K==1){
      om_inits <- sample.int(10, size = M, replace = TRUE)
      dim(init_vals$alpha) <- 1
      dim(init_vals$beta) <- 1
      dim(init_vals$w) <- 1
      init_vals$omega_k <- om_inits/sum(om_inits)
      dim(init_vals$omega_k) <- c(1,M)
    } 
    else{
      om_inits <- sample.int(10, size = M, replace = TRUE)
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
          //print(\"i = \", i);
          //print(\"Hs: \");
          //print(H_a[i]);
          //print(H_b[i]);
          //print(\"Hs_truncated: \");
          //print((H_a[i]*0.9943) + 0.00015);
          //print((H_b[i]*0.9943) + 0.00015);
          //print(alpha[k]);
          //print(beta[k]);
          log_weighted_mixt_cdf_a[k] += beta_proportion_lcdf((H_a[i] * 0.9943) + 0.00015| alpha[k], beta[k]);
          log_weighted_mixt_cdf_b[k] += beta_proportion_lcdf((H_b[i] * 0.9943) + 0.00015| alpha[k], beta[k]);
        }
        target += log_diff_exp(log_sum_exp(log_weighted_mixt_cdf_a), log_sum_exp(log_weighted_mixt_cdf_b));
      }
    }
    "
    w_vec_inits <- sample.int(100, size = K, replace = TRUE)
    init_vals <- list(
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
                            verbose=FALSE)
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