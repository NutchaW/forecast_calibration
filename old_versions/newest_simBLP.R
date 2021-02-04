##----------------------------individual bias correction----------------------------##

# individual bias correction by regressing Y on the mean of each point, 
# each Y is a realization of dis
# this is by definition of bias

# to get coefficient ( to biased correct values)
bc_slr<-function(Y, means){
  Yind<-Y[ind_est]
  means_ind<-means[ind_est, ]
  g1e<-lm(Yind~means_ind[,1])$coefficients
  g2e<-lm(Yind~means_ind[,2])$coefficients
  g3e<-lm(Yind~means_ind[,3])$coefficients
  coef_g<-as.matrix(rbind(g1e,g2e,g3e),ncol=2)
  return(coef_g)
}

##------------------------------Pooling Functions -----------------------------------##

wrapper_TLP <- function(Y, means, sds, weights, eval = TRUE){
  weights<-exp(weights)
  weights <- weights/sum(weights)
  ll <- log(dnorm(Y, mean = means, sd = sds) %*% weights)
  if(eval){
    return(sum(ll))
  } 
  else{
    return(ll)
  }
}

wrapper_bcBLP <- function(Y, means, sds, pars, coefs, eval = TRUE){
  ab     <- exp(pars[1:2])
  weights <- exp(pars[3:6])
  weights <- weights/sum(weights)
  
  bcY <- cbind(rep(1,nrow(Y)),Y) %*% coefs
  g      <- dnorm(bcY, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(bcY, mean = means, sd = sds) %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  ll     <- log(g) + log(g_beta)
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

wrapper_BLP <- function(Y, means, sds, pars, eval = TRUE){
  ab     <- exp(pars[1:2])
  weights <- exp(pars[-c(1,2)])
  weights <- weights/sum(weights)
  
  g      <- dnorm(Y, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds) %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  ll     <- log(g) + log(g_beta)
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

wrapper_bcBLP <- function(Y, means, sds, pars, coefs, eval = TRUE){
  ab     <- exp(pars[1:2])
  weights <- exp(pars[3:6])
  weights <- weights/sum(weights)
  
  bcY <- cbind(rep(1,nrow(Y)),Y) %*% coefs
  g      <- dnorm(bcY, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(bcY, mean = means, sd = sds) %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  ll     <- log(g) + log(g_beta)
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

# while function 
# ensemble
get_TLPpars<-function(means,sds){
  weights_TLP <- log(c(0.25, 0.25, 0.25,0.25))
#  while(TRUE){
    est <- optim(
      par     = weights_TLP,
      fn      = wrapper_TLP,
      Y       = Y[ind_est],
      means   = means[ind_est, ],
      sds   = sds[ind_est, ],
      method  = "Nelder-Mead",
      control = list(fnscale = -1, trace = TRUE)
    )
   weights_TLP <- exp(est$par)/sum(exp(est$par))
  return(weights_TLP)
}

get_BLPpars<-function(means,sds){
  weights_BLP <- log(c(0.25, 0.25, 0.25,0.25))
  ab          <- log(c(1,1))
  pars_BLP    <- c(ab, weights_BLP)
    est <- optim(
      par     = pars_BLP,
      fn      = wrapper_BLP,
      Y       = Y[ind_est],
      means   = means[ind_est, ],
      sds   = sds[ind_est, ],
      method  = "Nelder-Mead",
      control = list(fnscale = -1, trace = TRUE)
    )
    
    ab          <- exp(est$par[1:2])
    weights_BLP <- exp(est$par[-c(1,2)])/sum(exp(est$par[-c(1,2)]))
    pars_BLP    <- c(ab, weights_BLP)
  return(pars_BLP)
}

get_bcBLPpars<-function(means,sds){
  weights_bcBLP <- log(c(0.25, 0.25, 0.25,0.25))
  ab          <- log(c(1,1))
  pars_bcBLP    <- c(ab, weights_BLP)
  coefs<- bc_slr(Y[ind_est],means[ind_est,])
  est <- optim(
    par     = pars_bcBLP,
    fn      = wrapper_bcBLP,
    bcY       = Y[ind_est],
    means   = means[ind_est, ],
    sds   = sds[ind_est, ],
    coefs = coefs,
    method  = "Nelder-Mead",
    control = list(fnscale = -1, trace = TRUE)
  )
  ab          <- exp(est$par[1:2])
  weights_bcBLP <- exp(est$par[3:6])/sum(exp(est$par[3:6]))
  pars_bcBLP    <- c(ab, weights_bcBLP)
  return(list(pars_bcBLP,coefs))
}

##---------------------- PIT and log likelihood functions--------------------------------##
PIT_normal <- function(Y, mu, sd, var = FALSE, plot = FALSE,name){
  PIT <- pnorm(Y, mean = mu, sd = sd)
  
  if(var) return(var(PIT))
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 2), main = paste0("Density Forecast",name))
    abline( h = 1, lty = 2)
  }}


ll_normal         <- function(Y, mu, sd){
  ll <- log(dnorm(Y, mean = mu, sd = sd))
  return(ll)
}


wrapper_TLP_nonc <- function(Y, means, sds, weights, eval = TRUE, plot = FALSE, varPIT=FALSE){
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


wrapper_BLP_nonc  <- function(Y, means, sds, pars, eval = TRUE, plot = FALSE,varPIT=FALSE){
  
  ab     <- pars[1:2]
  weights <- pars[-c(1,2)]

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

wrapper_bcBLP <- function(Y, means, sds, pars, coefs, eval = TRUE,plot = FALSE,varPIT=FALSE){
  ab     <- exp(pars[1:2])
  weights <- exp(pars[3:6])
  weights <- weights/sum(weights)
  
  bcY <- cbind(rep(1,nrow(Y)),Y) %*% coefs
  g      <- dnorm(bcY, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(bcY, mean = means, sd = sds) %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  ll     <- log(g) + log(g_beta)
  PIT <- pbeta(pnorm(bcY, mean = means, sd = sds) %*% weights, shape1 = ab[1], shape2 = ab[2])
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "bcBLP")
    abline( h = 1, lty = 2)
  }  
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}