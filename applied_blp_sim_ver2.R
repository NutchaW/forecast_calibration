#---------------------------------------------- input functions ------------------------------------------------#
cdf_to_pdf <- function(cdf_vector){
  prob_vector <- c(cdf_vector[1],c(cdf_vector - dplyr::lag(cdf_vector))[-1])
  pdf <- prob_vector/sum(prob_vector)
  return(pdf)
}

#-------------------------------------------------- ensemble functions -----------------------------------------------#
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
  
  model_obj <- stan_model(model_code=model_type)
  pars <- rstan::optimizing(model_obj,data=input_data,seed=1234,hessian=T,init=initial_vals,
                            verbose=TRUE)
  #make a list of results
  params <- pars$par
  return(params)
}

## build ensemble pdfs
# make_pdfs <- function(ensemble_name, params, K,component_pdfs){
#   # generate some inputs
#   M <- unique(component_pdfs$model)
#   component_cdfs <- component_pdfs %>%
#     dplyr::group_by(model,index) %>%
#     dplyr::mutate(cdf_val=cumsum(binned_probs)) %>%
#     ungroup() %>%
#     dplyr::select(-"binned_probs")
#   comp_matrix <- component_pdfs %>%
#     tidyr::pivot_wider(id_cols =c(bin_start,bin_end,index), 
#                        names_from = model, 
#                        values_from = binned_probs) %>%
#     as.matrix()
#   comp_cdf_matrix <- component_cdfs %>%
#     dplyr::select("bin_start","bin_end","model","cdf_val","index") %>%
#     tidyr::pivot_wider(id_cols =c(bin_start,bin_end,index),
#                        names_from = model,
#                        values_from = binned_probs) %>%
#     as.matrix()
#   # frame for pdf
#   ensemble_pdf <- data.frame(comp_matrix[,1:4])
#   # build train and test cdf/pdf part
#   if(ensemble_name=="EW"|ensemble_name=="TLP"){
#     pdf_vals <- comp_matrix[,4:(4+M)] %*% params
#     ensemble_pdf$binned_probs <- pdf_vals
#     ensemble_pdf <- ensemble_pdf %>%
#       dplyr::group_by(index) %>%
#       dplyr::mutate(binned_probs=binned_probs/sum(binned_probs)) %>%
#       ungroup()
#   } 
#   else if(grepl("EW-", ensemble_name, fixed = TRUE)){
#     # ew-bmc or ew-blp
#     ## extract params
#     omega <- rep(1/M,M)
#     munu <- params[1:(K*2)]
#     w <- params[(1+(2*K)):length(params)]
#     # equally-weighed component models (ew ensemble)
#     H <- comp_cdf_matrix[,4:(4+M)] %*% omega
#     # build matrices for each k
#     pbeta_mat <- matrix(NA,ncol=K,nrow=nrow(comp_cdf_matrix))
#     for (k in 1:K){
#       assign(paste0("ab",k),c(munu[k]*munu[k+K],(1-munu[k])*munu[k+K]))
#       pbeta_mat[,k] <- pbeta( H, shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
#     }
#     # turn pbeta_mat into a combined mixture (pbeta_mat %*% w)
#     ensemble_pdf$cdf_vals <- pbeta_mat %*% w
#     # then turn  a combined mixture into pdf using the cdf to pdf function
#     ensemble_pdf <- ensemble_pdf %>%
#       dplyr::group_by(index) %>%
#       dplyr::mutate(binned_probs=cdf_to_pdf(cdf_vals)) %>%
#       ungroup() %>%
#       dplyr::select(-"cdf_vals")
#   } else {
#     # BMC or blp
#     ## extract params
#     munu <- params[1:(K*2)]
#     omega <- params[1+(K*2):(length(params)-(K+1))]
#     w <- params[(1+(2*K)+(M*K)):length(params)]
#     # build matrices for each k
#     pbeta_mat <- matrix(NA,ncol=K,nrow=nrow(comp_cdf_matrix))
#     for (k in 1:K){
#       assign(paste0("ab",k),c(munu[k]*munu[k+K],(1-munu[k])*munu[k+K]))
#       assign(paste0("omega",k),c())
#       for (m in 1:M){
#         assign(paste0("omega",k), c(get(paste0("omega",k)),omega[((m-1)*K)+k]))
#       }
#       H <- comp_cdf_matrix[,4:(4+M)] %*% c(get(paste0("omega",k)))
#       pbeta_mat[,k] <- pbeta( H, shape1 = get(paste0("ab",k))[1], shape2 = get(paste0("ab",k))[2])
#       }
#     # turn pbeta_mat into a combined mixture (pbeta_mat %*% w)
#     ensemble_pdf$cdf_vals <- pbeta_mat %*% w
#     # then turn  a combined mixture into pdf using the cdf to pdf function
#     ensemble_pdf <- ensemble_pdf %>%
#       dplyr::group_by(index) %>%
#       dplyr::mutate(binned_probs=cdf_to_pdf(cdf_vals)) %>%
#       ungroup() %>%
#       dplyr::select(-"cdf_vals")
#   }
#   return(ensemble_pdf)
# }

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
