## simulating Y
n <- 1e5
coefs <- c(a0 = 1, a1 = 1, a2 = 1, a3 = 1.1)
vars  <- matrix(rnorm(n = n*4), nrow = n, dimnames = list(NULL, c("x0", "x1", "x2", "x3")))
Y       <- vars %*% coefs + rnorm(n = n)

# Scenario 0: Simulation setup 
means0 <- sds0 <- matrix(NA, nrow = length(Y), ncol = 3)

for(i in 1:ncol(means0)){
    ind <- 1:4 %in% c(1, i+1)
    means0[,i]  <- vars[, ind] %*% coefs[ind]
    sds0[,i]    <- sqrt(1 + sum(coefs[!ind]^2))
}

## summary statistics
pits <- matrix(
    c(
        pnorm(Y, means0[,1], sds0[,1]),
        pnorm(Y, means0[,2], sds0[,2]),
        pnorm(Y, means0[,3], sds0[,3])
    ),
    ncol=3, byrow=FALSE
)
hist(pits[,1])
hist(pits[,2])
hist(pits[,3])

bias <- means0-cbind(Y, Y, Y)
colMeans(bias)
apply(bias, 2, FUN=var)

## Scenario 2: model 2 is underconfident, all unbiased
means2 <- sds2 <- matrix(NA, nrow = length(Y), ncol = 3)

for(i in 1:ncol(means2)){
    ind <- 1:4 %in% c(1, i+1)
    means2[,i]  <- (vars[, ind]) %*% coefs[ind]
    sds2[,i]    <- sqrt(1 + sum(coefs[!ind]^2))
}

## adding mis-calibration? kind of like shrinking to the middle?
#means2[,2]<- means2[,2]+sample(c(1,-1), n, replace=TRUE)
sds2[,1]<- sds2[,1]-0.1 ## subtracting variance to make it overconfident

## these two lines above have similar effects, both create overconfidence


## summary statistics
pits <- matrix(
    c(
        pnorm(Y, means2[,1], sds2[,1]),
        pnorm(Y, means2[,2], sds2[,2]),
        pnorm(Y, means2[,3], sds2[,3])
    ),
    ncol=3, byrow=FALSE
)
hist(pits[,1])
hist(pits[,2]) ## overconfident (dist'n too narrow): PIT is a smile
hist(pits[,3])

bias <- means2-cbind(Y, Y, Y)
colMeans(bias)
apply(bias, 2, FUN=var)
apply(means2, 2, FUN=var)

## Scenario 3: model 1 overconfident, all unbiased
means3 <- sds3 <- matrix(NA, nrow = length(Y), ncol = 3)

for(i in 1:ncol(means3)){
    ind <- 1:4 %in% c(1, i+1)
    means3[,i]  <- (vars[, ind]) %*% coefs[ind]
    sds3[,i]    <- sqrt(1 + sum(coefs[!ind]^2))
}

sds3[,1]<- sds3[,1]+0.5 ## adding variance to make it underconfident

## summary statistics
pits <- matrix(
    c(
        pnorm(Y, means3[,1], sds3[,1]),
        pnorm(Y, means3[,2], sds3[,2]),
        pnorm(Y, means3[,3], sds3[,3])
    ),
    ncol=3, byrow=FALSE
)
hist(pits[,1]) ## underconfident (too wide dist'n): PIT is a frown
hist(pits[,2]) 
hist(pits[,3])

bias <- means3-cbind(Y, Y, Y)
colMeans(bias)
apply(bias, 2, FUN=var)
