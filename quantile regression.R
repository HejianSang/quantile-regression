library(quantreg)
library(mvtnorm)
####################################################################################################
# The simulation studies were conducted to compare the performance of calibrated
# quantile estimator to those of direct estimator and difference estimator.
# Two finite populations of size N = 1000 were generated from bivariate
# normal distribution respectively. The correlation between two variables in
# the first population is 0.9 and in the second population is 0.6.



N=1000
pho1=0.6
pho2=0.9
mu=c(0,0)
sigma1=matrix(c(1,pho1,pho1,1),2,2)
sigma2=matrix(c(1,pho2,pho2,1),2,2)
pop1=rmvnorm(N,mean=mu,sigma=sigma1)
pop2=rmvnorm(N,mean=mu,sigma=sigma2)
