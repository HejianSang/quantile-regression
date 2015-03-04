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

samples=function(pop,n)
{
  N=dim(pop)[1]
  index=sample(1:N,n,replace=FALSE)
  S=pop
  S[!(1:N) %in% index,2]=NA
  S
}
A=samples(pop1,100)



#### This function in put a variable a and the points of segment that is a vector.
### This function return the index which group the a belongs
group=function(a,point)
{
  p=length(point)
  if(a<=point[1])
    return(1)
  for(i in 2:p)
  {
    if(a>point[i-1] & a<=point[i])
      return(i)
  }
  if(a>point[p])
    return(p+1)
}


category1=function(y,G)
{
  temp=NULL
  t=dim(y)[2]
  for(i in 1:t)
    temp=rbind(temp,quantile(y[,i],probs=seq(1/G,1-1/G,by=1/G),na.rm=T))
  ### seprate points
  yc=y
  for(i in 1:dim(y)[1])
    for(j in 1:dim(y)[2])
      if(!is.na(yc[i,j]))
      yc[i,j]=group(y[i,j],temp[j,])
  return(list(group=yc,point=temp))
}

### Function to find the poosible patterns in sample. 
### This help to reduce the heavy of computation
### Input the sample data which is the sample data we have grouped the miss place we use NA
### Return a complete matrix which is all possible patterns
possible_pat=function(sample_data)
{
  G=max(sample_data,na.rm=T)
  t=dim(sample_data)[2]
  pat=NULL
  for( i in 1:dim(sample_data)[1])
  {
    x=as.numeric(sample_data[i,])
    j=which(is.na(x))
    r=matrix(rep(x,G),ncol=t,byrow=T)
    r[,j]=1:G
    pat=rbind(pat,r)
  }
  pat
}

count=function(y_prime,rotate.cate)
{
  y_prime[is.na(y_prime)]=0
  temp=apply(rotate.cate,1,function(x)
    {
    x[is.na(x)]=0
    all(x==y_prime)
  })
  sum(temp)
}


EM=function(maxiter,joint_prob,rotate.cate,cate)
{
  G=max(rotate.cate,na.rm=T)
  N=dim(cate)[1]
  time=dim(cate)[2]
  p=dim(joint_prob)[1]
  weight=1
  for(t in 1:maxiter)
  {
    joint=sapply(1:p,function(i) 
    {
      y_prime=as.numeric(joint_prob[i,1:time])
      nt=count(y_prime,rotate.cate)*weight
      for(j in 1:time)
      {
        if(joint_prob[i,time+1]==0)
          break
        index=1:time
        y_prime=as.numeric(joint_prob[i,1:time])
        y_prime[j]=NA
        index=index[-j]
        if(count(y_prime,rotate.cate)!=0)
        {
          nplus=count(y_prime,rotate.cate)*weight*joint_prob[i,time+1]/ifelse(
            sum(joint_prob[sapply(joint_prob[,index],function(x) all(x==y_prime[index])),time+1])==0,1,
            sum(joint_prob[sapply(joint_prob[,index],function(x) all(x==y_prime[index])),time+1]))
          
        }
        else nplus=0
        nt=nt+nplus
      }
      return(nt)
    })
    #     if(sum((joint-joint_prob$prob)^2)<0.1)
    #       break
    joint_prob$prob=joint/sum(joint)
    cat("total=",sum(joint),"\n")
    cat("t=",t,"\n")
  }
  return(joint_prob)
}

em.est=function(tau,A,M,G)
{
  n=sum(!is.na(A[,2]))
  N=length(A[,1])
  X.pop=category1(A,10)$group[,1]
    index=sample(1:N,n,replace=FALSE)
    sample_data=A
    sample_data[!(1:N)%in% index,1]=NA
    cate.all=category1(sample_data,10)
    cate=cate.all$group
    rotate=sample_data[apply(is.na(sample_data),1,sum)<2,]
    rotate.cate=cate[apply(is.na(cate),1,sum)<2,]
    pat=data.frame(possible_pat(rotate.cate))
    joint_prob=data.frame(pat,prob=rep(1/dim(pat)[1],dim(pat)[1]))
    X.pop=sapply(A[,1],function(x) group(x,cate.all$point[1,]))
    margin.X=as.numeric(table(X.pop)/N)
    joint.prob=EM(M,joint_prob,rotate.cate,cate)
  #' Do calibration
  #' Use X.margin
  cali.weight=margin.X/(sapply(1:G,function(x) sum(joint.prob[joint.prob[,1]==x,"prob"])))
  margin.y=NULL
  for(j in 1:G)
  {
    m=0
    for(i in 1:G)
    {
      m=m+cali.weight[i]*sum(joint.prob[joint.prob[,1]==i&joint.prob[,2]==j,"prob"])
    }
    margin.y=c(margin.y,m)
  }
   
      
  Dist=cumsum(margin.y)
  index=max((1:G)[Dist<tau])
  left=tau-Dist[index]
  #'kernel smooth for group index +1
  d=A[cate[,2]==index+1,2]
  d=d[!is.na(A[cate[,2]==index+1,2])]
  temp=bkde(d,kernel="norm")
  diff=temp$x[2]-temp$x[1]
  min((temp$x[min((1:401)[temp$x>d[1]]):401])[cumsum(temp$y[min((1:401)[temp$x>d[1]]):401])*diff>left])
}

em.est(0.5,A,20,10)
