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

###############################################################
#' Here is the simple idea
#' We have SRS samples form Population Y
#' And we observe all X
#' First we category the sample y to be 10 groups
#' Then we use SRS to get same sample size for X
#'  Called sample x. We also category x to be 10 groups
#'  Then we use EM algorith to calculate the joint prob P(x=i,y=j)
#'  Calibration is we know that the marginal prob for X
#'  So sum P(x=i,y=j) for x is equal to p_i




