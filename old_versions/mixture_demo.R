library(pipeR)
library(Matrix)
library(tidyverse)
library(cowplot)
library(data.table)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(xtable)
library(here)
library(stats)
library(mixtools)
library(cowplot)
library(ggforce)
library(grid)
library(rmutil)
library(here)
library(knitr)
library(kableExtra)
library(rstan)
library(rjags)

set.seed(1234)

###### ----- try example ---------######
# mu <- c(-2.75, 2.75);
# sigma <- c(1, 1);
# lambda <- 0.4
# N <- 1000
# z <- rbinom(N, 1, lambda) + 1;
# y <- rnorm(N, mu[z], sigma[z]);
# # write file?
# rstan_options(auto_write = TRUE)
# stan_rdump(c("N", "y"), file="mix.data.R")
# gauss_mix <- "
# data {
#   int<lower = 0> N;
#   vector[N] y;
# }
# 
# parameters {
#   vector[2] mu;
#   real<lower=0> sigma[2];
#   real<lower=0, upper=1> theta;
# }
# 
# model {
#   sigma ~ normal(0, 2);
#   mu ~ normal(0, 2);
#   theta ~ beta(5, 5);
#   for (n in 1:N)
#     target += log_mix(theta,
#                       normal_lpdf(y[n] | mu[1], sigma[1]),
#                       normal_lpdf(y[n] | mu[2], sigma[2]));
# }
# "
# 
# input_data <- read_rdump("mix.data.R")
# degenerate_fit <- stan(model_code = gauss_mix, data=input_data,chains=4, seed=1234, refresh=2000)
# print(degenerate_fit)

############# --------------------------- simulation ---------------------------  ###############
# data

n <- 1e5
coefs <- c(a0 = 1, a1 = 1, a2 = 1, a3 = 1.1)
vars  <- matrix(rnorm(n = n*4), nrow = n, dimnames = list(NULL, c("x0", "x1", "x2", "x3")))

Y       <- vars %*% coefs + rnorm(n = n)

means0 <- sds0 <- matrix(NA, nrow = length(Y), ncol = 3)

for(i in 1:ncol(means0)){
  ind <- 1:4 %in% c(1, i+1)
  means0[,i]  <- vars[, ind] %*% coefs[ind]
  sds0[,i]    <- sqrt(1 + sum(coefs[!ind]^2))
}

# make density forecasts and component PITs
comp_forecasts <- dnorm(Y, means0, sds0)
comp_forecasts[comp_forecasts==0]<-0.00000001
comp_pits <- pnorm(Y, means0, sds0)
comp_pits[comp_pits==0]<-0.00000001
# other parameters 
K <- 2 # number of mixture components
M <- 3 # number of component models
alpha_w <- c(1.0,1.0)
alpha_omega <- c(1.0,1.0,1.0)
# write file?
input_data <- list(n=n, K=K, M=M, comp_forecasts=comp_forecasts, comp_pits=comp_pits,
                   alpha_w=alpha_w,alpha_omega=alpha_omega)
## code
sim_mix <- "
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

fit_s1 <- stan_model(model_code = sim_mix)
params <- optimizing(fit_s1,data=input_data,seed=1234,hessian=T)


# # ------------------ investigate mixture ------------------------- #
# 
# # data
# n <- 1e5
# p <- c(0.3,0.2,0.5)
# mu <- c(-2,0,3)
# sigma <- rep(1, 3)
# Y <- rnormmix(n, p, mu, sigma)
# comp_forecasts2 <- cbind(dnorm(Y, -2, 1),dnorm(Y, 0, 1),dnorm(Y, 3, 1))
# comp_pits2 <- cbind(pnorm(Y, -2, 1),pnorm(Y, 0, 1),pnorm(Y, 3, 1))
# 
# # other parameters 
# K <- 2 # number of mixture components
# M <- 3 # number of component models
# alpha_w <- c(1.0,1.0)
# alpha_omega <- c(1.0,1.0,1.0)
# # write file?
# input_data <- list(n=n, K=K, M=M, comp_forecasts=comp_forecasts, comp_pits=comp_pits,
#                    alpha_w=alpha_w,alpha_omega=alpha_omega)
# ## code
# sim_mix <- "
# data {
# int<lower=0> n;
# int<lower=0> K;
# int<lower=0> M;
# vector<lower=0>[K] alpha_w;
# vector<lower=0>[M] alpha_omega;
# matrix[n,M] comp_forecasts;
# matrix[n,M] comp_pits;
# }
# 
# transformed data {
# matrix[n,M] log_f = log(comp_forecasts);
# }
# 
# parameters {
# vector<lower=0, upper=1>[K] mu;
# vector<lower=0>[K] nu;
# simplex[M] omega_k[K];
# simplex[K] w;
# }
# 
# model {
# matrix[n, K] H;
# vector[K] log_weighted_mixt_density;
# vector[M] inner_log_mixt_density;
# vector[K] log_w = log(w);
# vector[M] log_omega_k[K];
# // priors
# for (k in 1:K) {
# mu[k] ~ beta(2, 2);
# nu[k] ~ gamma(0.1, 0.1);
# omega_k[k] ~ dirichlet(alpha_omega);
# log_omega_k[k] = log(omega_k[k]);
# H[, k] = comp_pits*omega_k[k];
# }
# w ~ dirichlet(alpha_w);
# for (i in 1:n) {
# log_weighted_mixt_density = log_w;
# for (k in 1:K) {
# inner_log_mixt_density = log_omega_k[k];
# for (m in 1:M) {
# inner_log_mixt_density[m] += log_f[i, m];
# }
# log_weighted_mixt_density[k] += log_sum_exp(inner_log_mixt_density);
# log_weighted_mixt_density[k] += beta_proportion_lpdf(H[i, k] | mu[k], nu[k]);
# }
# target += log_sum_exp(log_weighted_mixt_density);
# }
# }
# "
# 
# fit <- stan_model(model_code = sim_mix)
# params <- optimizing(fit_s1,data=input_data,seed=1234,hessian=T)
# 
