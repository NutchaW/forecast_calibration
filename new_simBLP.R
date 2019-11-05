# Scenario 1:
# - assume all models are unbiased.
# - apply TLP, SLP, component-specific SLP (cSLP), BLP, cBLP 
#Functions -------------------------------------------------------------------
wrapper_TLP <- function(Y, means, sds, weights, eval = TRUE, 
                        plot = FALSE, varPIT = FALSE){
  if(any(weights > 1) || any(weights < 0)) return(-1e8)
  weights <- weights/sum(weights)
  ll <- log(dnorm(Y, mean = means, sd = sds) %*% weights)
  PIT <- pnorm(Y, mean = means, sd = sds) %*% weights
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "TLP")
    abline( h = 1, lty = 2)
  }
  if(varPIT) return(var(PIT))
  if(eval){
    return(sum(ll))
  } 
  else{
    return(ll)
  }
}

wrapper_SLP <- function(Y, means, sds, pars, eval = TRUE, 
                        plot = FALSE, varPIT = FALSE){
  
  con     <- pars[1]
  weights <- pars[-1]
  
  if(con <= 0) return(-1e8)
  if(any(weights > 1) || any(weights < 0)) return(-1e8)
  
  weights <- weights/sum(weights)
  
  means_SLP <- (matrix(rep(Y, ncol(means)), ncol = ncol(means)) - means)/con
  
  ll <- log(1/con * dnorm(means_SLP, mean = 0, sd = sds) %*% weights)
  
  PIT <- pnorm(means_SLP, mean = 0, sd = sds) %*% weights
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "SLP")
    abline( h = 1, lty = 2)
  }  
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

wrapper_BLP <- function(Y, means, sds, pars, eval = TRUE, 
                        plot = FALSE, varPIT = FALSE){
  
  ab     <- pars[1:2]
  weights <- pars[-c(1,2)]
  
  if(any(ab <= 0)) return(-1e8)
  if(any(weights > 1) || any(weights < 0)) return(-1e8)
  
  weights <- weights/sum(weights)
  
  g      <- dnorm(Y, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds) %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(pnorm(Y, mean = means, sd = sds) %*% weights, shape1 = ab[1], shape2 = ab[2])
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "BLP")
    abline( h = 1, lty = 2)
  }  
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

wrapper_NBLP <- function(Y, means, sds, pars, eval = TRUE, 
                        plot = FALSE, varPIT = FALSE){
  
  abc     <- pars[1:3]
  weights <- pars[-c(1:3)]
  
  if(any(abc[1:2] <= 0)) return(-1e8)
  if(any(abc[3] < 0)) return(-1e8)
  if(any(weights > 1) || any(weights < 0)) return(-1e8)
  
  weights <- weights/sum(weights)
  
  g      <- dnorm(Y, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds) %*% weights, 
                  shape1 = abc[1], shape2 = abc[2],ncp=abc[3])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(pnorm(Y, mean = means, sd = sds) %*% weights, shape1 = abc[1], shape2 = abc[2],ncp=abc[3])
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "Non-central BLP")
    abline( h = 1, lty = 2)
  }  
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

wrapper_indSLP <- function(Y, means, sds, pars, eval = TRUE, 
                        plot = FALSE, varPIT = FALSE){
  con     <- pars[1:3]
  weights <- pars[-(1:3)]
  if(any(con <= 0)) return(-1e8)
  if(any(weights > 1) || any(weights < 0)) return(-1e8)
  weights <- weights/sum(weights)
  means_scaled <- (matrix(rep(Y, ncol(means)), ncol = ncol(means)) - means)
  means_SLP <- matrix(c(means_scaled[,1]/con[1],means_scaled[,2]/con[2],means_scaled[,3]/con[3]),
                      ncol = ncol(means))
  #ll <- log(1/con * dnorm(means_SLP, mean = 0, sd = sds) %*% weights)
  ll <- log(cbind(dnorm(means_SLP[,1], mean = 0, sd = sds[,1])/con[1],
              dnorm(means_SLP[,2], mean = 0, sd = sds[,2])/con[2],
              dnorm(means_SLP[,3], mean = 0, sd = sds[,3])/con[3]) %*% weights)

  PIT <- pnorm(means_SLP, mean = 0, sd = sds) %*% weights
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "cSLP")
    abline( h = 1, lty = 2)
  }  
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

wrapper_indBLP <- function(Y, means, sds, pars, eval = TRUE){
  
  ab     <- pars[1:2]

  if(any(ab <= 0)) return(-1e8)
  g      <- dnorm(Y, mean = means, sd = sds) 
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds), 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
  if(eval){
    return(sum(ll))
  }else{
    return(ll)
  }
}


wrapper_BBLP <- function(Y, means, sds, pars, cab, eval = TRUE, 
                         plot = FALSE, varPIT = FALSE){
  comp_bpars<-cab
  # normal probs
  comp_nprobs<-cbind(pnorm(Y, mean = means[,1], sd = sds[,1]),
                     pnorm(Y, mean = means[,2], sd = sds[,2]),
                     pnorm(Y, mean = means[,3], sd = sds[,3]))
  # get beta probs
  comp_bprobs<-cbind(pbeta(pnorm(Y, mean = means[,1], sd = sds[,1]), 
                                  shape1 = comp_bpars[1,1], shape2 = comp_bpars[1,2]),
                            pbeta(pnorm(Y, mean = means[,2], sd = sds[,2]), 
                                  shape1 = comp_bpars[2,1], shape2 = comp_bpars[2,2]),
                            pbeta(pnorm(Y, mean = means[,3], sd = sds[,3]), 
                                  shape1 = comp_bpars[3,1], shape2 = comp_bpars[3,2]))
  ab     <- pars[1:2]
  weights <- pars[-c(1,2)]
  
  if(any(ab <= 0)) return(-1e8)
  if(any(weights > 1) || any(weights < 0)) return(-1e8)
  
  weights <- weights/sum(weights)
  f1 <- as.matrix(dnorm(Y, mean = means, sd = sds))
  f2<- as.matrix(dbeta(comp_nprobs,shape1 = comp_bpars[,1], shape2 = comp_bpars[,2]))
  g      <- (f1*f2) %*% weights
  g_beta <- dbeta(comp_bprobs %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(comp_bprobs %*% weights, shape1 = ab[1], shape2 = ab[2])
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "cBLP")
    abline( h = 1, lty = 2)
  }  
  
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}


# while function 

# Estimate weights for TLP
# change back
get_TLPpars<-function(means,sds){
  weights_TLP <- c(0.3, 0.3, 0.4)
  while(TRUE){
    
    est <- optim(
      par     = weights_TLP,
      fn      = wrapper_TLP,
      Y       = Y[ind_est],
      means   = means[ind_est, ],
      sds   = sds[ind_est, ],
      method  = "BFGS",
      control = list(fnscale = -1, trace = TRUE)
    )
    
    weights_TLP <- est$par/sum(est$par)
    if(abs(sum(est$par) - 1) < 0.01) break
  }
  return(weights_TLP)
}

# Estimate weights and c (= con) for SLP

get_SLPpars<-function(means,sds){
  weights_SLP <- c(0.3, 0.3, 0.4)
  con         <- 1
  pars_SLP    <- c(con, weights_SLP)
  while(TRUE){
    
    est <- optim(
      par     = pars_SLP,
      fn      = wrapper_SLP,
      Y       = Y[ind_est],
      means   = means[ind_est, ],
      sds   = sds[ind_est, ],
      method  = "BFGS",
      control = list(fnscale = -1, trace = TRUE)
    )
    
    con         <- est$par[1]
    weights_SLP <- est$par[-1]/sum(est$par[-1])
    pars_SLP    <- c(con, weights_SLP)
    if(abs(sum(est$par[-1]) - 1) < 0.01) break
  }
  return(pars_SLP)
}

get_indSLPpars<-function(means,sds){
  weights_SLP <- c(0.3, 0.3, 0.4)
  con         <- c(1,1,1)
  pars_SLP    <- c(con, weights_SLP)
  
  # Estimate weights and c (= con) for SLP
  while(TRUE){
    
    est <- optim(
      par     = pars_SLP,
      fn      = wrapper_indSLP,
      Y       = Y[ind_est],
      means   = means[ind_est, ],
      sds   = sds[ind_est, ],
      method  = "BFGS",
      control = list(fnscale = -1, trace = TRUE)
    )
    
    con         <- est$par[1:3]
    weights_SLP <- est$par[-c(1:3)]/sum(est$par[-c(1:3)])
    pars_SLP    <- c(con, weights_SLP)
    if(abs(sum(est$par[-c(1:3)]) - 1) < 0.01) break
  }
  return(pars_SLP)
}

get_BLPpars<-function(means,sds){
  weights_BLP <- c(0.3, 0.3, 0.4)
  ab          <- c(1,1)
  pars_BLP    <- c(ab, weights_BLP)
  
  # Estimate weights for BLP
  while(TRUE){
    
    est <- optim(
      par     = pars_BLP,
      fn      = wrapper_BLP,
      Y       = Y[ind_est],
      means   = means[ind_est, ],
      sds   = sds[ind_est, ],
      method  = "BFGS",
      control = list(fnscale = -1, trace = TRUE)
    )
    
    ab          <- est$par[1:2]
    weights_BLP <- est$par[-c(1,2)]/sum(est$par[-c(1,2)])
    pars_BLP    <- c(ab, weights_BLP)
    if(abs(sum(est$par[-c(1,2)]) - 1) < 0.01) break
  }
  return(pars_BLP)
}

# get_NBLPpars<-function(means,sds){
#   weights_nBLP <- c(0.3, 0.3, 0.4)
#   abc          <- c(1,1,0)
#   pars_nBLP    <- c(abc, weights_nBLP)
#   
#   # Estimate weights for BLP
#   while(TRUE){
#     
#     est <- optim(
#       par     = pars_nBLP,
#       fn      = wrapper_NBLP,
#       Y       = Y[ind_est],
#       means   = means[ind_est, ],
#       sds   = sds[ind_est, ],
#       method  = "BFGS",
#       # control = list(fnscale = -1, trace = TRUE)
#       control = list(fnscale = -1,trace = TRUE)
#       
#     )
#     
#     abc          <- est$par[1:3]
#     weights_nBLP <- est$par[-c(1:3)]/sum(est$par[-c(1:3)])
#     pars_BLP    <- c(abc, weights_nBLP)
#     if(est$par[c(1:2)]>0||est$par[3]>=0) break
#   }
#   return(pars_nBLP)
# }


get_indBLPpars<-function(means,sds,modelnum){
  ab          <- c(1,1)
  while(TRUE){
    est <- optim(
      par     = ab,
      fn      = wrapper_indBLP,
      Y       = Y[ind_est],
      means   = means[ind_est,modelnum],
      sds   = sds[ind_est,modelnum],
      method  = "BFGS",
      control = list(fnscale = -1, trace = TRUE)
    )
    ab          <- est$par
    #if(sum(est$par[-1]) < Inf) break # check if results the same
    if(est$par>0||sum(est$par) < Inf) break
  }
  return(ab)
}

get_BBLPpars<-function(means,sds,cab){
  weights_BBLP <- c(0.3, 0.3, 0.4)
  ab          <- c(1,1)
  pars_BBLP    <- c(ab, weights_BBLP)
  while(TRUE){
    
    est <- optim(
      par     = pars_BBLP,
      fn      = wrapper_BBLP,
      Y       = Y[ind_est],
      means   = means[ind_est, ],
      sds   = sds[ind_est, ],
      cab   = cab,
      method  = "BFGS",
      control = list(fnscale = -1, trace = TRUE)
    )
    ab          <- est$par[1:2]
    weights_BBLP <- est$par[-c(1,2)]/sum(est$par[-c(1,2)])
    pars_BBLP    <- c(ab, weights_BBLP)
    if(abs(sum(est$par[-c(1,2)]) - 1) < 0.01) break
  }
  return(pars_BBLP)
}


# PIT values of a normal distribution
PIT_normal <- function(Y, mu, sd, var = FALSE, plot = FALSE,name){
  PIT <- pnorm(Y, mean = mu, sd = sd)
  
  if(var) return(var(PIT))
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 2), main = paste0("Density Forecast",name))
    abline( h = 1, lty = 2)
  }}
# PIT values of a beta distribution
PIT_beta <- function(Y, mu, sd, ab, var = FALSE, plot = FALSE,name){
  PIT <- pbeta(pnorm(Y, mean = mu, sd = sd), shape1 = ab[1], shape2 = ab[2])
  if(var) return(var(PIT))
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = paste0("Beta Density Forecast",name))
    abline( h = 1, lty = 2)
  }  
}
# single bin log score wrapper
ll_normal         <- function(Y, mu, sd){
  ll <- log(dnorm(Y, mean = mu, sd = sd))
  return(ll)
}

ll_beta         <- function(Y,mu,sd, ab){
  probs<-pbeta(pnorm(Y, mean = mu, sd = sd), shape1 = ab[1], shape2 = ab[2])
  ll <- log(dnorm(Y, mean = mu, sd = sd))+log(dbeta(probs, shape1 = ab[1], shape2 = ab[2]))
  return(ll)
}
