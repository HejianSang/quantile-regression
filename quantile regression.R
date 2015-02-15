library(quantreg)
library(mvtnorm)
#'The simulation studies were conducted to compare the performance of calibrated
#' quantile estimator to those of direct estimator and difference estimator.
#' Two finite populations of size N = 1000 were generated from bivariate
#' normal distribution respectively. The correlation between two variables in
#' the first population is 0.9 and in the second population is 0.6.



N=1000
pho1=0.6
pho2=0.9
mu=c(0,0)
sigma1=matrix(c(1,pho1,pho1,1),2,2)
sigma2=matrix(c(1,pho2,pho2,1),2,2)
pop1=data.frame(rmvnorm(N,mean=mu,sigma=sigma1))
names(pop1)=c("X","Y")

pop2=data.frame(rmvnorm(N,mean=mu,sigma=sigma2))
names(pop2)=c("X","Y")

#'To use simple random sampling to get n=100 sample
n=100


#' Direct Method
Direct=function(pop,n,tau)
{
  N=dim(pop)[1]
  rho=cor(pop)[1,2]
  D=sapply(1:1000,function(o) {
    A=pop[sample(1:N,n,replace=FALSE),]
    weights=rep(N/n,n)
    r.regression=rq(Y~X-1,tau=tau,data=A,weights=weights)
    y.hat=r.regression$coef*A[,1]
    y=sapply(A[,1],function(x) qnorm(tau,mean=rho*x,sd=sqrt(1-rho^2)))
    sum((y-y.hat)^2)/n
  })
  
  mean(D)
}

tau=c(0.1,0.3,0.5,0.7,0.9)
mse.d1=sapply(tau,function(tau) Direct(pop1,n,tau))
mse.d2=sapply(tau,function(tau) Direct(pop2,n,tau))







