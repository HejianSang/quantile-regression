library(quantreg)
library(mvtnorm)
library(plyr)
library(ggplot2)
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

#' Use weights to get the quantile for tau
#' y is a vector of saples
#' tau is the quantile probability
#' weight is the vector for each sample response 

weighted.quantile=function(y,tau,weight)
{
  data=cbind(y=y,weight=weight)
  data=data[order(data[,1]),]
  data=data.frame(data)
  total.weight=sum(weight)
  f.hat=cumsum(data$weight)/total.weight
  indice=which(f.hat>=tau)
  return(data$y[min(indice)])
}

#' Direct Method and
#' We wan to use difference estimator 
#' Difference estimator also use the auxilary information\
#' So we want to compare with callibartion method
compare=function(pop,n,tau,tau0)
{
  N=dim(pop)[1]
  rho=cor(pop)[1,2]
  pop.regression=rq(Y~X-1,tau=tau0,data=pop)
  beta.star=pop.regression$coef
  pop.q.mean=mean(pop[,1]*beta.star)
  D=sapply(1:1000,function(o) {
    A=pop[sample(1:N,n,replace=FALSE),]
    weights=rep(N/n,n)
    r.regression=rq(Y~X-1,tau=tau0,data=A,weights=weights)
    q.direct=weighted.quantile(A[,2],tau,weights)
    q.true=quantile(pop[,2],probs=tau)
    mse1=(q.direct-q.true)^2
    q.diff=q.direct+(pop.q.mean-mean(beta.star*A[,1]))
    mse2=(q.diff-q.true)^2
    B=sum(r.regression$coef*pop[,1]*(pop[,1]*beta.star-pop.q.mean))/sum((pop[,1]*beta.star-pop.q.mean)^2)
    q.emp=q.direct+(pop.q.mean-mean(beta.star*A[,1]))*B
    mse3=(q.emp-q.true)^2
    return(c(mse1,mse2,mse3))
  })
  
  apply(D,1,mean)
}

tau=c(0.1,0.3,0.5,0.7,0.9)
mse=sapply(tau,function(tau) compare(pop1,n,tau,0.3))
mse[1,]/mse[3,]
mse[2,]/mse[3,]

mse1=sapply(tau,function(tau) compare(pop1,n,tau,0.7))
mse1[1,]/mse1[3,]
mse1[2,]/mse1[3,]
#########################################################
mse=sapply(tau,function(tau) compare(pop2,n,tau,0.3))
mse[1,]/mse[3,]
mse[2,]/mse[3,]

mse1=sapply(tau,function(tau) compare(pop2,n,tau,0.7))
mse1[1,]/mse1[3,]
mse1[2,]/mse1[3,]

#' We want to consider the different weights d_i
#' we want to use q(x_i,beta_hat) to replace q(x_i,beta_prob)
#' because they are asympototic equal.
compare.v2=function(pop,weights,tau,tau0)
{
  N=dim(pop)[1]
  rho=cor(pop)[1,2]
  n=length(weights)
#   pop.regression=rq(Y~X-1,tau=tau0,data=pop)
#   beta.star=pop.regression$coef
#   pop.q.mean=mean(pop[,1]*beta.star)
  D=sapply(1:1000,function(o) {
    A=pop[sample(1:N,n,replace=FALSE),]
    weights=weights
    r.regression=rq(Y~X,tau=tau0,data=A,weights=weights)
    q.direct=weighted.quantile(A[,2],tau,weights)
    q.true=quantile(pop[,2],probs=tau)
    mse1=(q.direct-q.true)^2
    q.diff=q.direct+(sum(r.regression$coef*pop[,1])-sum(weights*(r.regression$coef*A[,1])))/N
    mse2=(q.diff-q.true)^2
#     B=sum(r.regression$coef*pop[,1]*(pop[,1]*beta.star-pop.q.mean))/sum((pop[,1]*beta.star-pop.q.mean)^2)
#     q.emp=q.direct+(pop.q.mean-mean(beta.star*A[,1]))*B
    mse3=mse2
    return(c(mse1,mse2,mse3))
  })
  
  apply(D,1,mean)
}

tau=c(0.1,0.3,0.5,0.7,0.9)
mse=sapply(tau,function(tau) compare.v2(pop1,rep(N/n,n),tau,0.3))
mse[1,]/mse[3,]
############################################################################
#' To estimate the population quantile at tau1, 
#' what tau we use to do callibrartion is most effective?
tau0=seq(0,1,0.05)
tau=0.5
mse.check=sapply(tau0,function(tau0) compare(pop1,n,tau,tau0))
mse.cb=data.frame(tau0=tau0,mse.c=mse.check[3,])
png(filename="tau_5.png")
qplot(tau0,mse.c,data=mse.cb,geom="point")
dev.off()
tau=0.7
mse.check=sapply(tau0,function(tau0) compare(pop1,n,tau,tau0))
mse.cb=data.frame(tau0=tau0,mse.c=mse.check[3,])
png(filename="tau_7.png")
qplot(tau0,mse.c,data=mse.cb,geom="point")
dev.off()

#########################################################################
#' Variance estimation
#' We have two method to approach the variance 
#' The first one is nonparametric method







