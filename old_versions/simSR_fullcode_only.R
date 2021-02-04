library("pipeR");library("dplyr");library("cdcfluview")
library(fGarch);library(Matrix)
library(quantreg);library(cowplot)
library(data.table);library(reshape2);library(gridExtra)
library(rmutil);library(ranger)
library(xgboost);library(splines);library(ggplot2)
library(xtable);library(here)
library(stats);library(mixtools);library(cowplot)
library(ggforce);library(grid);library(FluSight)
library(rmutil);library(here);library(stats)
library(mixtools);library(knitr);library(kableExtra)

source("./new_simBLP.R")
set.seed(123)
# -----------Scenario 0: Simulation setup -----------------------#
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

forecast_names <- c("f1", "f2", "f3", "TLP", "SLP", "BLP")

ind_est <- c(rep(TRUE, n/2), rep(FALSE, n/2))

# TLP (Traditional linear pool) combination -----------------------------------
# Simple linear combination of forecasts
weights_TLP0<-get_TLPpars(means0,sds0)
round(weights_TLP0, 3)

# SLP (Spread-adjusted linear pool) combination -------------------------------
par_SLP0<-get_SLPpars(means0,sds0)
weights_SLP0<-par_SLP0[-1]
con0<-par_SLP0[1]
round(weights_SLP0, 3)
round(con0, 3)

# component-specific SLP (Spread-adjusted linear pool) combination -------------------------------
par_indSLP0<-get_indSLPpars(means0,sds0)
weights_indSLP0<-par_SLP0[-c(1:3)]
ind_con0<-par_indSLP0[c(1:3)]
round(weights_indSLP0, 3)
round(ind_con0, 3)

# BLP (Beta-transformed linear pool) combination ------------------------------
par_BLP0<-get_BLPpars(means0,sds0)
weights_BLP0<-par_BLP0[-c(1:2)]
ab0<-par_BLP0[c(1:2)]
round(weights_BLP0, 3)
round(ab0, 3)

# NBLP (Beta-transformed linear pool) combination ------------------------------
par_NBLP0<-get_NBLPpars(means0,sds0)
weights_NBLP0<-par_NBLP0[-c(1:3)]
abc0<-par_NBLP0[c(1:3)]
round(weights_NBLP0, 3)
round(abc0, 3)

# component-specific BLP (Beta-transformed linear pool) combination ------------------------------
comp_bpars0<-matrix(NA,ncol=2,nrow=3)
for(i in 1:3){
  ab0<-get_indBLPpars(means0,sds0,i)
  comp_bpars0[i,1]<-ab0[1]
  comp_bpars0[i,2]<-ab0[2]
}
round(comp_bpars0, 3)

# BLP components in BLP
par_CBLP0<-get_CBLPpars(means0,sds0,comp_bpars0)
weights_CBLP0<-par_CBLP0[-c(1:2)]
abb0<-par_CBLP0[c(1:2)]
round(weights_CBLP0, 3)
round(abb0, 3)

# Evaluation ------------------------------------------------------------------

# Parameter values
table0 <- matrix(NA, nrow = 8, ncol = 9, 
                 dimnames = list(c(forecast_names[4:6],"cSLP","beta-f1","beta-f2","beta-f3","cBLP"), 
                                 c("w1", "w2", "w3", "c", "alpha", "beta",
                                   "c1","c2","c3")))
table0[1,] <- c(weights_TLP0, rep(NA, 6))
table0[2,] <- c(par_SLP0[-1], par_SLP0[1], rep(NA, 5))
table0[3,] <- c(par_BLP0[-c(1,2)], NA, par_BLP0[1:2],rep(NA, 3))
table0[4,] <- c(par_indSLP0[-c(1:3)],rep(NA, 3),par_indSLP0[1:3])
table0[5,] <- c(rep(NA, 4),comp_bpars0[1,c(1:2)],rep(NA, 3))
table0[6,] <- c(rep(NA, 4),comp_bpars0[2,c(1:2)],rep(NA, 3))
table0[7,] <- c(rep(NA, 4),comp_bpars0[3,c(1:2)],rep(NA, 3))
table0[8,] <- c(par_CBLP0[-c(1,2)], NA, par_CBLP0[1:2],rep(NA, 3))
table0<-round(table0, 2)

# var(PIT) 
vPIT0 <- matrix(NA,nrow=11,ncol = 1)
rownames(vPIT0) <- c(forecast_names,"cSLP","beta-f1","beta-f2","beta-f3","cBLP")
colnames(vPIT0) <- c("Variance")

for(i in 1:3) vPIT0[i] <- PIT_normal(Y = Y[!ind_est], mu = means0[!ind_est ,i], 
                                     sd = sds0[!ind_est ,i], var = TRUE)
vPIT0[4] <- wrapper_TLP(Y = Y[!ind_est], means = means0[!ind_est, ],
                        sds = sds0[!ind_est, ], weights = weights_TLP0, varPIT = TRUE)
vPIT0[5] <- wrapper_SLP(Y = Y[!ind_est], means = means0[!ind_est, ],
                        sds = sds0[ind_est, ], pars = par_SLP0, varPIT = TRUE)
vPIT0[6] <- wrapper_BLP(Y = Y[!ind_est], means = means0[!ind_est, ],
                        sds = sds0[!ind_est, ], pars = par_BLP0, varPIT = TRUE)
vPIT0[7] <- wrapper_indSLP(Y = Y[!ind_est], means = means0[!ind_est, ],
                           sds = sds0[!ind_est, ], pars = par_indSLP0, varPIT = TRUE)
for(i in 1:3) vPIT0[7+i] <- PIT_beta(Y = Y[!ind_est], mu = means0[!ind_est ,i], 
                                     sd = sds0[!ind_est ,i],ab=comp_bpars0[i,c(1:2)], var = TRUE)
vPIT0[11] <- wrapper_CBLP(Y = Y[!ind_est], means = means0[!ind_est, ],
                          sds = sds0[!ind_est, ], pars = par_CBLP0, cab=comp_bpars0, varPIT = TRUE)

# # Mean log score
# # # LS(f,y) = -log(f(y))
mls0 <- matrix(NA, nrow = 11, ncol = 2,
               dimnames = list(c(forecast_names,
                                 "cSLP","beta-f1","beta-f2","beta-f3","cBLP"),
                               c("Training", "Test")))
for(j in 1:2){
  ind <- if(j == 1){ind_est}else{!ind_est}
  for(i in 1:3) mls0[i, j]  <- mean(ll_normal(Y = Y[ind,], mu = means0[ind, i],
                                              sd = sds0[ind, i]))
  mls0[4, j] <- mean(wrapper_TLP(Y = Y[ind], means = means0[ind,],
                                 sds = sds0[ind,], weights = weights_TLP0, eval = FALSE))
  mls0[5, j] <- mean(wrapper_SLP(Y = Y[ind], means = means0[ind,],
                                 sds = sds0[ind,], pars = par_SLP0, eval = FALSE))
  mls0[6, j] <- mean(wrapper_BLP(Y = Y[ind], means = means0[ind,],
                                 sds = sds0[ind,], pars = par_BLP0, eval = FALSE))
  mls0[7, j] <- mean(wrapper_indSLP(Y = Y[ind], means = means0[ind,],
                                    sds = sds0[ind,], pars = par_indSLP0, eval = FALSE))
  for(i in 1:3) mls0[7+i, j]  <- mean(ll_beta(Y = Y[ind,], mu = means0[ind, i],
                                              sd = sds0[ind, i], ab=comp_bpars0[i,]))
  mls0[11, j] <- mean(wrapper_CBLP(Y = Y[ind], means = means0[ind,], sds = sds0[ind,],
                                   pars = par_CBLP0,cab=comp_bpars0, eval = FALSE))
}

# PITplot(0)
par(mfrow = c(2,3))
for(i in 1:3) PIT_normal(Y = Y[!ind_est], mu = means0[!ind_est ,i],
                         sd = sds0[!ind_est ,i], plot = TRUE, name = i)
wrapper_TLP(Y = Y[!ind_est], means = means0[!ind_est, ],
            sds = sds0[!ind_est, ], weights = weights_TLP0, plot = TRUE)
wrapper_SLP(Y = Y[!ind_est], means = means0[!ind_est, ],
            sds = sds0[!ind_est, ], pars = par_SLP0, plot = TRUE)
wrapper_BLP(Y = Y[!ind_est], means = means0[!ind_est, ],
            sds = sds0[!ind_est, ], pars = par_BLP0, plot = TRUE)
wrapper_indSLP(Y = Y[!ind_est], means = means0[!ind_est, ],
               sds = sds0[!ind_est, ], pars = par_indSLP0, plot = TRUE)
for(i in 1:3) PIT_beta(Y = Y[!ind_est], mu = means0[!ind_est ,i],
                       sd = sds0[!ind_est ,i],ab=comp_bpars0[i,c(1:2)], plot = TRUE, name = i)
wrapper_CBLP(Y = Y[!ind_est], means = means0[!ind_est, ],
             sds = sds0[!ind_est, ], pars = par_CBLP0, cab=comp_bpars0, plot = TRUE)

#------------ Scenario 2: some models are overconfident, others underconfident? (edited)---------- 
#
#
# Simulation setup ------------------------------------------------------------

coefs2 <- c(a0 = 1, a1 = 1, a2 = 1, a3 = 1.1)
means2 <- sds2 <- matrix(NA, nrow = length(Y), ncol = 3)

for(i in 1:ncol(means2)){
  ind <- 1:4 %in% c(1, i+1)
  means2[,i]  <- (vars[, ind]) %*% coefs2[ind]
  sds2[,i]    <- sqrt(1 + sum(coefs2[!ind]^2))
}

means2[,2]<- means2[,2]+sample(c(1,-1), n, replace=TRUE)

# TLP (Traditional linear pool) combination -----------------------------------
# Simple linear combination of forecasts
weights_TLP2<-get_TLPpars(means2,sds2)
round(weights_TLP2, 3)

# SLP (Spread-adjusted linear pool) combination -------------------------------
par_SLP2<-get_SLPpars(means2,sds2)
weights_SLP2<-par_SLP2[-1]
con2<-par_SLP2[1]
round(weights_SLP2, 3)
round(con2, 3)

# component-specific SLP (Spread-adjusted linear pool) combination -------------------------------
par_indSLP2<-get_indSLPpars(means2,sds2)
weights_indSLP2<-par_SLP2[-c(1:3)]
ind_con2<-par_indSLP2[c(1:3)]
round(weights_indSLP2, 3)
round(ind_con2, 3)

# BLP (Beta-transformed linear pool) combination ------------------------------
par_BLP2<-get_BLPpars(means2,sds2)
weights_BLP2<-par_BLP2[-c(1:2)]
ab2<-par_BLP2[c(1:2)]
round(weights_BLP2, 3)
round(ab2, 3)

# NBLP (Beta-transformed linear pool) combination ------------------------------
par_NBLP2<-get_NBLPpars(means2,sds2)
weights_NBLP2<-par_NBLP2[-c(1:3)]
abc2<-par_NBLP2[c(1:3)]
round(weights_NBLP2, 3)
round(abc2, 3)

# component-specific BLP (Beta-transformed linear pool) combination ------------------------------
comp_bpars2<-matrix(NA,ncol=2,nrow=3)
for(i in 1:3){
  ab2<-get_indBLPpars(means2,sds2,i)
  comp_bpars2[i,1]<-ab2[1]
  comp_bpars2[i,2]<-ab2[2]
}
round(comp_bpars2, 3)

# BLP components in BLP
par_CBLP2<-get_CBLPpars(means2,sds2,comp_bpars2)
weights_CBLP2<-par_CBLP2[-c(1:2)]
abb2<-par_CBLP2[c(1:2)]
round(weights_CBLP2, 3)
round(abb2, 3)

# Evaluation ------------------------------------------------------------------

# Parameter values
table2 <- matrix(NA, nrow = 8, ncol = 9, 
                 dimnames = list(c(forecast_names[4:6],"cSLP","beta-f1","beta-f2","beta-f3","cBLP"), 
                                 c("w1", "w2", "w3", "c", "alpha", "beta",
                                   "c1","c2","c3")))
table2[1,] <- c(weights_TLP2, rep(NA, 6))
table2[2,] <- c(par_SLP2[-1], par_SLP2[1], rep(NA, 5))
table2[3,] <- c(par_BLP2[-c(1,2)], NA, par_BLP2[1:2],rep(NA, 3))
table2[4,] <- c(par_indSLP2[-c(1:3)],rep(NA, 3),par_indSLP2[1:3])
table2[5,] <- c(rep(NA, 4),comp_bpars2[1,c(1:2)],rep(NA, 3))
table2[6,] <- c(rep(NA, 4),comp_bpars2[2,c(1:2)],rep(NA, 3))
table2[7,] <- c(rep(NA, 4),comp_bpars2[3,c(1:2)],rep(NA, 3))
table2[8,] <- c(par_CBLP2[-c(1,2)], NA, par_CBLP2[1:2],rep(NA, 3))
table2<-round(table2, 2)

# var(PIT) 

vPIT2 <- matrix(NA,nrow=11,ncol = 1)
rownames(vPIT2) <- c(forecast_names,"cSLP","beta-f1","beta-f2","beta-f3","cBLP")
colnames(vPIT2) <- c("Variance")
for(i in 1:3) vPIT2[i] <- PIT_normal(Y = Y[!ind_est], mu = means2[!ind_est ,i], 
                                     sd = sds2[!ind_est ,i], var = TRUE)
vPIT2[4] <- wrapper_TLP(Y = Y[!ind_est], means = means2[!ind_est, ],
                        sds = sds2[!ind_est, ], weights = weights_TLP2, varPIT = TRUE)
vPIT2[5] <- wrapper_SLP(Y = Y[!ind_est], means = means2[!ind_est, ],
                        sds = sds2[ind_est, ], pars = par_SLP2, varPIT = TRUE)
vPIT2[6] <- wrapper_BLP(Y = Y[!ind_est], means = means2[!ind_est, ],
                        sds = sds2[!ind_est, ], pars = par_BLP2, varPIT = TRUE)
vPIT2[7] <- wrapper_indSLP(Y = Y[!ind_est], means = means2[!ind_est, ],
                           sds = sds2[!ind_est, ], pars = par_indSLP2, varPIT = TRUE)
for(i in 1:3) vPIT2[7+i] <- PIT_beta(Y = Y[!ind_est], mu = means2[!ind_est ,i], 
                                     sd = sds2[!ind_est ,i],ab=comp_bpars2[i,c(1:2)], var = TRUE)
vPIT2[11] <- wrapper_CBLP(Y = Y[!ind_est], means = means2[!ind_est, ],
                          sds = sds2[!ind_est, ], pars = par_CBLP2, cab=comp_bpars2, varPIT = TRUE)

# # Mean log score
# # # LS(f,y) = -log(f(y))
mls2 <- matrix(NA, nrow = 11, ncol = 2,
               dimnames = list(c(forecast_names,
                                 "cSLP","beta-f1","beta-f2","beta-f3","cBLP"),
                               c("Training", "Test")))
for(j in 1:2){
  ind <- if(j == 1){ind_est}else{!ind_est}
  for(i in 1:3) mls2[i, j]  <- mean(ll_normal(Y = Y[ind,], mu = means2[ind, i],
                                              sd = sds2[ind, i]))
  mls2[4, j] <- mean(wrapper_TLP(Y = Y[ind], means = means2[ind,],
                                 sds = sds2[ind,], weights = weights_TLP2, eval = FALSE))
  mls2[5, j] <- mean(wrapper_SLP(Y = Y[ind], means = means2[ind,],
                                 sds = sds2[ind,], pars = par_SLP2, eval = FALSE))
  mls2[6, j] <- mean(wrapper_BLP(Y = Y[ind], means = means2[ind,],
                                 sds = sds2[ind,], pars = par_BLP2, eval = FALSE))
  mls2[7, j] <- mean(wrapper_indSLP(Y = Y[ind], means = means2[ind,],
                                    sds = sds2[ind,], pars = par_indSLP2, eval = FALSE))
  for(i in 1:3) mls2[7+i, j]  <- mean(ll_beta(Y = Y[ind,], mu = means2[ind, i],
                                              sd = sds2[ind, i], ab=comp_bpars2[i,]))
  mls2[11, j] <- mean(wrapper_CBLP(Y = Y[ind], means = means2[ind,], sds = sds2[ind,],
                                   pars = par_CBLP2,cab=comp_bpars2, eval = FALSE))
}

# PIT Histogramme
# 
par(mfrow = c(2,3))
for(i in 1:3) PIT_normal(Y = Y[!ind_est], mu = means2[!ind_est ,i],
                         sd = sds2[!ind_est ,i], plot = TRUE, name = i)
wrapper_TLP(Y = Y[!ind_est], means = means2[!ind_est, ],
            sds = sds2[!ind_est, ], weights = weights_TLP2, plot = TRUE)
wrapper_SLP(Y = Y[!ind_est], means = means2[!ind_est, ],
            sds = sds2[!ind_est, ], pars = par_SLP2, plot = TRUE)
wrapper_BLP(Y = Y[!ind_est], means = means2[!ind_est, ],
            sds = sds2[!ind_est, ], pars = par_BLP2, plot = TRUE)
wrapper_indSLP(Y = Y[!ind_est], means = means2[!ind_est, ],
               sds = sds2[!ind_est, ], pars = par_indSLP2, plot = TRUE)
for(i in 1:3) PIT_beta(Y = Y[!ind_est], mu = means2[!ind_est ,i],
                       sd = sds2[!ind_est ,i],ab=comp_bpars2[i,c(1:2)], plot = TRUE, name = i)
wrapper_CBLP(Y = Y[!ind_est], means = means2[!ind_est, ],
             sds = sds2[!ind_est, ], pars = par_CBLP2, cab=comp_bpars2, plot = TRUE)

# Scenario 3: 
coefs3 <- c(a0 = 1, a1 = 1, a2 = 1, a3 = 1.1)
means3 <- sds3 <- matrix(NA, nrow = length(Y), ncol = 3)

for(i in 1:ncol(means3)){
  ind <- 1:4 %in% c(1, i+1)
  means3[,i]  <- (vars[, ind]) %*% coefs3[ind]
  sds3[,i]    <- sqrt(1 + sum(coefs3[!ind]^2))
}

sds3[,1]<- sds3[,1]+0.5

# TLP (Traditional linear pool) combination -----------------------------------
# Simple linear combination of forecasts
weights_TLP3<-get_TLPpars(means3,sds3)
round(weights_TLP3, 3)

# SLP (Spread-adjusted linear pool) combination -------------------------------
par_SLP3<-get_SLPpars(means3,sds3)
weights_SLP3<-par_SLP3[-1]
con3<-par_SLP3[1]
round(weights_SLP3, 3)
round(con3, 3)

# component-specific SLP (Spread-adjusted linear pool) combination -------------------------------
par_indSLP3<-get_indSLPpars(means3,sds3)
weights_indSLP3<-par_SLP3[-c(1:3)]
ind_con3<-par_indSLP3[c(1:3)]
round(weights_indSLP3, 3)
round(ind_con3, 3)

# BLP (Beta-transformed linear pool) combination ------------------------------
par_BLP3<-get_BLPpars(means3,sds3)
weights_BLP3<-par_BLP3[-c(1:2)]
ab3<-par_BLP3[c(1:2)]
round(weights_BLP3, 3)
round(ab3, 3)


# component-specific BLP (Beta-transformed linear pool) combination ------------------------------
comp_bpars3<-matrix(NA,ncol=2,nrow=3)
for(i in 1:3){
  ab3<-get_indBLPpars(means3,sds3,i)
  comp_bpars3[i,1]<-ab3[1]
  comp_bpars3[i,2]<-ab3[2]
}
round(comp_bpars3, 3)


# BLP components in BLP
par_CBLP3<-get_CBLPpars(means3,sds3,comp_bpars3)
weights_CBLP3<-par_CBLP3[-c(1:2)]
abb3<-par_CBLP3[c(1:2)]
round(weights_CBLP3, 3)
round(abb3, 3)

# Evaluation ------------------------------------------------------------------

# Parameter values
table3 <- matrix(NA, nrow = 8, ncol = 9, 
                 dimnames = list(c(forecast_names[4:6],"cSLP","beta-f1","beta-f2","beta-f3","cBLP"), 
                                 c("w1", "w2", "w3", "c", "alpha", "beta",
                                   "c1","c2","c3")))
table3[1,] <- c(weights_TLP3, rep(NA, 6))
table3[2,] <- c(par_SLP3[-1], par_SLP3[1], rep(NA, 5))
table3[3,] <- c(par_BLP3[-c(1,2)], NA, par_BLP3[1:2],rep(NA, 3))
table3[4,] <- c(par_indSLP3[-c(1:3)],rep(NA, 3),par_indSLP3[1:3])
table3[5,] <- c(rep(NA, 4),comp_bpars3[1,c(1:2)],rep(NA, 3))
table3[6,] <- c(rep(NA, 4),comp_bpars3[2,c(1:2)],rep(NA, 3))
table3[7,] <- c(rep(NA, 4),comp_bpars3[3,c(1:2)],rep(NA, 3))
table3[8,] <- c(par_CBLP3[-c(1,2)], NA, par_CBLP3[1:2],rep(NA, 3))
table3<-round(table3, 2)

# var(PIT) 
vPIT3 <- matrix(NA,nrow=11,ncol = 1)
rownames(vPIT3) <- c(forecast_names,"cSLP","beta-f1","beta-f2","beta-f3","cBLP")
colnames(vPIT3) <- c("Variance")
for(i in 1:3) vPIT3[i] <- PIT_normal(Y = Y[!ind_est], mu = means3[!ind_est ,i], 
                                     sd = sds3[!ind_est ,i], var = TRUE)
vPIT3[4] <- wrapper_TLP(Y = Y[!ind_est], means = means3[!ind_est, ],
                        sds = sds3[!ind_est, ], weights = weights_TLP3, varPIT = TRUE)
vPIT3[5] <- wrapper_SLP(Y = Y[!ind_est], means = means3[!ind_est, ],
                        sds = sds3[ind_est, ], pars = par_SLP3, varPIT = TRUE)
vPIT3[6] <- wrapper_BLP(Y = Y[!ind_est], means = means3[!ind_est, ],
                        sds = sds3[!ind_est, ], pars = par_BLP3, varPIT = TRUE)
vPIT3[7] <- wrapper_indSLP(Y = Y[!ind_est], means = means3[!ind_est, ],
                           sds = sds3[!ind_est, ], pars = par_indSLP3, varPIT = TRUE)
for(i in 1:3) vPIT3[7+i] <- PIT_beta(Y = Y[!ind_est], mu = means3[!ind_est ,i], 
                                     sd = sds3[!ind_est ,i],ab=comp_bpars3[i,c(1:2)], var = TRUE)
vPIT3[11] <- wrapper_CBLP(Y = Y[!ind_est], means = means3[!ind_est, ],
                          sds = sds3[!ind_est, ], pars = par_CBLP3, cab=comp_bpars3, varPIT = TRUE)
vPIT3<-round(vPIT3, 3)
# # Mean log score
# # # LS(f,y) = -log(f(y))
mls3 <- matrix(NA, nrow = 11, ncol = 2,
               dimnames = list(c(forecast_names,
                                 "cSLP","beta-f1","beta-f2","beta-f3","cBLP"),
                               c("Training", "Test")))
for(j in 1:2){
  ind <- if(j == 1){ind_est}else{!ind_est}
  for(i in 1:3) mls3[i, j]  <- mean(ll_normal(Y = Y[ind,], mu = means3[ind, i],
                                              sd = sds3[ind, i]))
  mls3[4, j] <- mean(wrapper_TLP(Y = Y[ind], means = means3[ind,],
                                 sds = sds3[ind,], weights = weights_TLP3, eval = FALSE))
  mls3[5, j] <- mean(wrapper_SLP(Y = Y[ind], means = means3[ind,],
                                 sds = sds3[ind,], pars = par_SLP3, eval = FALSE))
  mls3[6, j] <- mean(wrapper_BLP(Y = Y[ind], means = means3[ind,],
                                 sds = sds3[ind,], pars = par_BLP3, eval = FALSE))
  mls3[7, j] <- mean(wrapper_indSLP(Y = Y[ind], means = means3[ind,],
                                    sds = sds3[ind,], pars = par_indSLP3, eval = FALSE))
  for(i in 1:3) mls3[7+i, j]  <- mean(ll_beta(Y = Y[ind,], mu = means3[ind, i],
                                              sd = sds3[ind, i], ab=comp_bpars3[i,]))
  mls3[11, j] <- mean(wrapper_CBLP(Y = Y[ind], means = means3[ind,], sds = sds3[ind,],
                                   pars = par_CBLP3,cab=comp_bpars3, eval = FALSE))
}

# PIT Histogramme
par(mfrow = c(2,3))
for(i in 1:3) PIT_normal(Y = Y[!ind_est], mu = means3[!ind_est ,i],
                         sd = sds3[!ind_est ,i], plot = TRUE, name = i)
wrapper_TLP(Y = Y[!ind_est], means = means3[!ind_est, ],
            sds = sds3[!ind_est, ], weights = weights_TLP3, plot = TRUE)
wrapper_SLP(Y = Y[!ind_est], means = means3[!ind_est, ],
            sds = sds3[!ind_est, ], pars = par_SLP3, plot = TRUE)
wrapper_BLP(Y = Y[!ind_est], means = means3[!ind_est, ],
            sds = sds3[!ind_est, ], pars = par_BLP3, plot = TRUE)
wrapper_indSLP(Y = Y[!ind_est], means = means3[!ind_est, ],
               sds = sds3[!ind_est, ], pars = par_indSLP3, plot = TRUE)
for(i in 1:3) PIT_beta(Y = Y[!ind_est], mu = means3[!ind_est ,i],
                       sd = sds3[!ind_est ,i],ab=comp_bpars3[i,c(1:2)], plot = TRUE, name = i)
wrapper_CBLP(Y = Y[!ind_est], means = means3[!ind_est, ],
             sds = sds3[!ind_est, ], pars = par_CBLP3, cab=comp_bpars3, plot = TRUE)