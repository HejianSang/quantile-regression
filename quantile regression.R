library(quantreg)
library(mvtnorm)
library(plyr)
library(ggplot2)
library(KernSmooth)
library(rootSolve)
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
    q.diff=q.direct+(mean(r.regression$coef*pop[,1])-mean(r.regression$coef*A[,1]))
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
  n=length(weights)
#   pop.regression=rq(Y~X-1,tau=tau0,data=pop)
#   beta.star=pop.regression$coef
#   pop.q.mean=mean(pop[,1]*beta.star)
  D=sapply(1:1000,function(o) {
    A=pop[sample(1:N,n,replace=FALSE),]
    weights=weights
    r.regression=rq(Y~X-1,tau=tau0,data=A,weights=weights)
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
mse=sapply(tau,function(tau) compare.v2(pop1,rep(N/n,n),tau,tau))
mse[1,]/mse[3,]

########################################################################
#' We also want to use more pre quantiles as calibrartion like tao0=c(0.1,0.2,...,0.9)
#' We can also try the moment calibrartion
#' We compare them with the regression estimators or calibrartion estimators
compare.v3=function(pop,tau)
{
  N=dim(pop)[1]
  D=sapply(1:1000,function(o) {
    A=pop[sample(1:N,n,replace=FALSE),]
    weights=rep(N/n,n)
    coef=NULL
    tau0=seq(0.1,0.9,0.2)
    for(i in 1:length(tau0))
      coef=cbind(coef,rq(Y~X,tau=tau0,data=A,weights=weights)$coef)
    q.direct=weighted.quantile(A[,2],tau,weights)
    q.true=quantile(pop[,2],probs=tau)
    mse1=(q.direct-q.true)^2
    q.diff=q.direct+(sum(coef[1,3]+coef[2,3]*pop[,1])-sum(weights*(coef[1,3]+coef[2,3]*A[,1])))/N
    mse2=(q.diff-q.true)^2
#     q.more=q.direct
#     for( i in 1:length(tau0))
#       q.more=q.more+(sum(coef[1,i]+coef[2,i]*pop[,1])-sum(weights*(coef[1,i]+coef[2,i]*A[,1])))/N
    q.ratio=weighted.quantile(pop[,1],tau,rep(1,N))*weighted.quantile(
      A[,2],tau,weights)/weighted.quantile(A[,1],tau,weights)
    mse3=(q.ratio-q.true)^2
    X=cbind(rep(1,N),pop[,1],pop[,1]^2,pop[,1]^3,pop[,1]^4)
    x=cbind(rep(1,n),A[,1],A[,1]^2,A[,1]^3,A[,1]^4)
    x=data.frame(x)
    data=cbind(x,y=A[,2])
    beta=lm(y~X2+X3+X4+X5,data=data)$coef
    q.mom=q.direct+sum((apply(X,2,mean)-apply(x,2,mean))*beta)
    mse4=(q.mom-q.true)^2
    return(c(mse1,mse2,mse3,mse4))
  })
  
  apply(D,1,mean)
}  
  
tau=c(0.1,0.3,0.5,0.7,0.9)
mse3=sapply(tau,function(tau) compare.v3(pop1,tau))
mse3[3,]/mse3[1,]
mse3[2,]/mse3[1,] 
mse3[4,]/mse3[1,] 
 
  
  
  
  







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



###################################################################
#' function to use callibrartion to estimate the quantile
Fw=function(y,theta,weight)
{
  sum(weight[y<=as.numeric(theta)])/sum(weight)
  
}
#####################################################
#' max emperical likelihood to get the weights
max.weight=function(q,d,q.N.bar)
{
  n=length(q)
  f=function(lambda)
  {
    c(sum(d/(lambda[1]+lambda[2]*q))+1,sum((d/(lambda[1]+lambda[2]*q)*q))+q.N.bar)
  }
  lambda=multiroot(f,start=c(-1,-1),maxiter=100)$root
  -d/(lambda[1]+q*lambda[2])
#   w=-d/(lambda[1]+q*lambda[2])
#   w[w<0]=0
#   w
}
#########################################################################
#' Variance estimation
#' We have two method to approach the variance 
#' The first one is nonparametric method
#' We need use kernel to estimate the density at the quantiel we estimate
#' We prefer the normal kernel 
#' of course we can try other kernel

var.kernel=function(pop,n,tau,tau0)
{
  N=dim(pop)[1]
  #' SRS samples
  A=pop[sample(1:N,n,replace=FALSE),]
  #' quantile regression
  weights=rep(N/n,n)
  PI=rep(n/N,n)
  r.regression=rq(Y~X,tau=tau0,data=A,weights=weights)
  q.N.bar=mean(r.regression$coef[1]+r.regression$coef[2]*pop[,1])
  q.direct=weighted.quantile(A[,2],tau,weights)
  theta.w=q.direct+(sum(r.regression$coef[1]+r.regression$coef[2]*pop[,1])-sum(weights*(r.regression$coef[1]+r.regression$coef[2]*A[,1])))/N
  q=r.regression$coef[1]+r.regression$coef[2]*A[,1]
  Den=matrix(0,2,2)
  Num=matrix(0,2,1)
  for(i in  1:n)
  {
    Den=Den+weights[i]*matrix(c(1,q[i],q[i],q[i]^2),2,2)
    Num=Num+weights[i]*matrix(c(1,q[i]),2,1)*(A[i,2]<theta.w)
  }
  C=solve(Den)%*%Num
  z=(A[,2]<=theta.w)-matrix(c(rep(1,n),q),n,2)%*%C
  PIJ=matrix(n*(n-1)/N/(N-1),n,n)
  diag(PIJ)=1/weights
  V.hat=var(z)*(1-n/N)/n
#   for(i in 1:n)
#     for(j in 1:n)
#       V.hat=V.hat+(PIJ[i,j]-PI[i]*PI[j])*z[i]*z[j]/(PIJ[i,j]*PI[i]*PI[j])
  #' to estimate the density of y
#   density=bkde(A[,2],kernel="normal")
#   index=min(which(density$x>theta.w))
#   f.theta=(density$y[index]+density$y[index-1])/2
  temp=matrix(c(rep(1,N),pop[,1]),N,2)
  q.N.bar=mean(temp%*%t(t( r.regression$coef)))
  h=2*sqrt(V.hat)
  w=max.weight(q,weights,q.N.bar)
  f.theta1=(Fw(A[,2],theta.w+h,w)-Fw(A[,2],theta.w-h,w))/(2*h)
#' f.theta could be 0
#' That will cause the problem

  density=bkde(A[,2],kernel="normal")
  index=min(which(density$x>theta.w))
  f.theta=(density$y[index]+density$y[index-1])/2
  V.theta=V.hat/f.theta^2
  return(list(theta.w=theta.w,V.hat=V.hat,V.theta=V.theta,f.theta1=f.theta1))
  
}

var.kernel(pop1,100,0.5,0.5)

#' To calculate the bias of the variance estimate
#' We repeat the procedure 5000 times
#' We get th 5000 theta estimations and variance estimations
#' we use the mean of variance estimations as expected variance
#' we use 50000 theta to calculate v(theta)
r1=sapply(1:5000,function(o) var.kernel(pop1,100,0.5,0.5))
E.V.hat1=mean(unlist(r1[3,]))
V.true1=var(unlist(r1[1,]))
(E.V.hat1/V.true1-1)*100



######################################################
#' Wooddruff method
#' Fw inverse quantile function with weights
Fw.quant=function(w,tau,tau0,A,pop)
{
  N=dim(pop)[1]
  r.regression=rq(Y~X,tau=tau0,data=A,weights=w)
  q.direct=weighted.quantile(A[,2],tau,w)
  theta.w=q.direct+(sum(r.regression$coef[1]+r.regression$coef[2]*pop[,1])-sum(w*(r.regression$coef[1]+r.regression$coef[2]*A[,1])))/N
  theta.w
}







wooddruff=function(pop,n,tau,tau0)
{
  N=dim(pop)[1]
  #' SRS samples
  A=pop[sample(1:N,n,replace=FALSE),]
  #' quantile regression
  weights=rep(N/n,n)
  PI=rep(n/N,n)
  r.regression=rq(Y~X,tau=tau0,data=A,weights=weights)
  q.direct=weighted.quantile(A[,2],tau,weights)
  theta.w=q.direct+(sum(r.regression$coef[1]+r.regression$coef[2]*pop[,1])-sum(weights*(r.regression$coef[1]+r.regression$coef[2]*A[,1])))/N
  q=r.regression$coef[1]+r.regression$coef[2]*A[,1]
  Den=matrix(0,2,2)
  Num=matrix(0,2,1)
  for(i in  1:n)
  {
    Den=Den+weights[i]*matrix(c(1,q[i],q[i],q[i]^2),2,2)
    Num=Num+weights[i]*matrix(c(1,q[i]),2,1)*(A[i,2]<theta.w)
  }
  C=solve(Den)%*%Num
  z=(A[,2]<theta.w)-matrix(c(rep(1,n),q),n,2)%*%C
  PIJ=matrix(n*(n-1)/N/(N-1),n,n)
  diag(PIJ)=1/weights
  V.hat=var(z)*(1-n/N)/n
  tau.L=as.numeric(tau-1.96*sqrt(V.hat))
  tau.U=as.numeric(tau+1.96*sqrt(V.hat))
  q.direct1=weighted.quantile(A[,2],tau.L,weights)
  theta.L=q.direct1+(sum(r.regression$coef[1]+r.regression$coef[2]*pop[,1])-sum(w*(r.regression$coef[1]+r.regression$coef[2]*A[,1])))/N
  q.direct2=weighted.quantile(A[,2],tau.U,weights)
  theta.U=q.direct2+(sum(r.regression$coef[1]+r.regression$coef[2]*pop[,1])-sum(w*(r.regression$coef[1]+r.regression$coef[2]*A[,1])))/N
  V.theta=(theta.U-theta.L)^2/4
  return(list(theta.w=theta.w,V.hat=V.hat,V.theta=V.theta))
}

wooddruff(pop1,100,0.5,0.5)

r=sapply(1:5000,function(o) wooddruff(pop2,100,0.2,0.2))
E.V.hat=mean(unlist(r[3,]))
V.true=var(unlist(r[1,]))
(E.V.hat/V.true-1)*100
############################################################