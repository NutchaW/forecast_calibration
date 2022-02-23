a <- 6.875
b <- 2066.612
omega <- c(0.98,0.02)
x<- comp_forecasts3 %*% omega
cdf5<- pbeta(x,shape1=a,shape2=b)
pdf5<- dbeta(x,shape1=a,shape2=b)

plot(x,cdf5[,1])
plot(x,pdf5[,1])
