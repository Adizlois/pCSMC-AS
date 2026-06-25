# Author: Alfonso Diz-Lois Palomares
# Email: adpaloma@uio.no
source("Aux_functions.R")
library(matrixStats)

###IMPLEMENTATION OF MARGINAL PCSMC WITH ANCESTOR SAMPLING
pCSMC_marginal = function(y,xstar,R_t,B,
                 hosp_prob,AS,alpha=NULL,beta=NULL)
{
  nT = dim(xstar)[1]
  p = dim(xstar)[2]-1
  if (is.null(alpha)) alpha <- 2
  if (is.null(beta))  beta  <- 4
  #Theta prior (gamma distribution)
  alphas=rep(alpha,p)
  betas=rep(beta,p)
  #
  xsim = array(0,c(nT,p+1,B)) #Number of infections. Each "p" represents how many of the infections of that 
  #day come from each of the previous days. Last dimension is just the sum
  suffs=array(0,c(2,nT,p,B)) # 2 sufficients: counts and exposure (R_t*I_{t-m})
  
  xsim[,,1] = xstar #Set reference as first particle
  w = rep(0,B)
  #Seeding
  xsim[1,1:p,-1] = 2
  xsim[1,p+1,-1]= 2*p
  
  w = w+ dbinom(y[1],xsim[1,p+1,],hosp_prob,log=TRUE)
  w[is.na(w)]<-0
  ############
  #Update sufficient
  suffs[1,,,]<-alpha
  suffs[2,,,]<-beta

  for(s in 2:nT)
  {

    a1 = 1
    #Ancestor sampling step
    if ((s>2)&(AS)){
    wtilde = w
    xstar_array <- array(
      rep(xstar[s:nT,,drop=F], times = B),
      dim = c(dim(xstar[s:nT,,drop=F]), B)
    )
    xsim_all <- abind::abind(xsim[1:(s-1),,,drop=F], xstar_array, along = 1)
    kmax=min(p,s-1)
    for (j in 1:kmax){
      wtilde=wtilde+dnbinom(xstar[s,j],size=suffs[1,s-1,j,],prob=suffs[2,s-1,j,]/(suffs[2,s-1,j,]+R_t[s] * xsim_all[s-j, p+1,]),log=T)
      }
    
    if (s<nT){
      suffs_as=suffs[,s-1,,]
    for (tt in (s+1):nT){
      kmax=min(p,tt-1)
      for (j in 1:kmax){
        suffs_as[1,j,]=suffs_as[1,j,]+xsim_all[tt-1,j,]
        if (tt-1-j>0){
        suffs_as[2,j,]=suffs_as[2,j,]+R_t[tt-1]*xsim_all[(tt-1-j),p+1,]
        }
        wtilde=wtilde+dnbinom(xsim_all[tt,j,],size=suffs_as[1,j,],prob=suffs_as[2,j,]/(suffs_as[2,j,]+R_t[tt]*xsim_all[(tt-j),p+1,]),log=T)
      }
     }
    }
    finite <- is.finite(wtilde)
    m <- max(wtilde[finite])
    wtilde[finite] <- exp(wtilde[finite] - m)
    wtilde[!finite]<-0
    a1 = sample(1:B,1,prob=wtilde)
    } ##End of Ancestor Sampling
    
    w = exp(w-max(w))
    #Resampling
    a = c(a1,sample(1:B,B-1,prob=w,replace=TRUE))
    xsim[1:(s-1),,] = xsim[1:(s-1),,a] 
    suffs[,1:(s-1),,]=suffs[,1:(s-1),,a]
    w = rep(0,B)
    #Forward model
    
    kmax <- min(p, s-1)
    #Update sufficient 2 at time s
    suffs[2,s,1:kmax,]=suffs[2,s-1,1:kmax,]+R_t[s]*xsim[(s-1):(s-kmax),p+1,]
    for(j in 1:kmax){
      xsim[s,j,-1 ] <- rnbinom(B-1,size=suffs[1,s-1,j,-1],prob=suffs[2,s-1,j,-1]/(suffs[2,s-1,j,-1]+R_t[s]*xsim[(s-j),p+1,-1]))
      }
    
    xsim[s:nT,,1]<-xstar[s:nT,]
    
    
    xsim[s, ,-1][is.na(xsim[s, ,-1])] <- 0
    xsim[s,p+1,-1]=apply(xsim[s,1:p,-1],2,sum,na.rm=T)
    #Update sufficient 1 at time s
    suffs[1,s,1:kmax,]=suffs[1,s-1,1:kmax,]+xsim[s,1:kmax,]
    w = w+ dbinom(y[s],xsim[s,p+1,],hosp_prob,log=TRUE)
    w[is.infinite(w)]<--100 #Arbritrary
  }
  #Sample reference at last time point
  w = exp(w-max(w))
  k = sample(1:B,1,prob=w)
  return(list(latent=xsim[,,k],param=samptheta1(xsim =xsim[,,k],R_t=R_t)))
}


##SMC used in the first iteration to initialize the Xs without Ancestor Sampling.
SMC_marginal = function(y,R_t,B,p,
                          hosp_prob,alpha=NULL,beta=NULL)
{
  nT =length(y)
  if (is.null(alpha)) alpha <- 2
  if (is.null(beta))  beta  <- 4
  #Theta prior (gamma distribution)
  alphas=rep(alpha,p)
  betas=rep(beta,p)
  #
  xsim = array(0,c(nT,p+1,B)) #Number of infections. Each "p" represents how many of the infections of that 
  #day come from each of the previous days. Last dimension is just the sum
  suffs=array(0,c(2,nT,p,B)) # 2 sufficients: counts and exposure (R_t*I_{t-m})
  w = rep(0,B)
  #Seeding
  xsim[1,1:p,] = 2
  xsim[1,p+1,]= 2*p
  
  w = w+ dbinom(y[1],xsim[1,p+1,],hosp_prob,log=TRUE)
  w[is.na(w)]<-0
  #For the first p points, we will just use the prior two options: based on the reference (as in Finke%Thiery) or pure prior
  ############
  #Update sufficient
  suffs[1,,,]<-alpha
  suffs[2,,,]<-beta
  
  for(s in 2:nT)
  {

    w = exp(w-max(w))
    
    ##Resampling
    a = c(sample(1:B,B,prob=w,replace=TRUE))
    xsim[1:(s-1),,] = xsim[1:(s-1),,a] 
    suffs[,1:(s-1),,]=suffs[,1:(s-1),,a]
    w = rep(0,B)
    #Forward model
    
    kmax <- min(p, s-1)
    suffs[2,s,1:kmax,]=suffs[2,s-1,1:kmax,]+R_t[s]*xsim[(s-1):(s-kmax),p+1,]
    for(j in 1:kmax){
      xsim[s,j,] <- rnbinom(B,size=suffs[1,s-1,j,],prob=suffs[2,s-1,j,]/(suffs[2,s-1,j,]+R_t[s]*xsim[(s-j),p+1,]))
    }
    xsim[s,p+1,]=apply(xsim[s,1:p,],2,sum,na.rm=T)
    suffs[1,s,1:kmax,]=suffs[1,s-1,1:kmax,]+xsim[s,1:kmax,]
    w = w+ dbinom(y[s],xsim[s,p+1,],hosp_prob,log=TRUE)
    w[is.infinite(w)]<--100 #Arbritrary
  }
  w = exp(w-max(w))
  k = sample(1:B,1,prob=w)
  return(list(latent=xsim[,,k],param=samptheta1(xsim =xsim[,,k],R_t=R_t)))
}

##MAIN FUNCTION
pSMCGibbs_marginal = function(y,R_t,hosp_prob,p,B=100,M=100,sseed=10)
{
  set.seed(sseed)
  nT = length(y)
  xsimM = array(NA,c(M+1,nT,p+1))
  thetaM = array(NA,c(M+1,p))
  # Initialize xsim with SMC
  res = SMC_marginal(y=y,R_t=R_t,B=B,p=p,
                       hosp_prob=hosp_prob)
  xsim<-res$latent
  xsimM[1,,]<-xsim 
  for(m in 1:M)
  {
    cat("\r",m, " ")
    #Sample x using CSMC
      res = pCSMC_marginal(y,xstar=xsim,R_t,B,
                           hosp_prob,AS=T)
    xsim=res$latent
    xsim[is.na(xsim)]<-0
    xsimM[m+1,,] = xsim
    
    thetaM[m+1,]=res$param
  }
  list(xsim=xsimM,thetaM=thetaM)
}
###########EXAMPLE USAGE

#The script assumes the synthetic data is available as syn_data.RDS. 
#Otherwise, run the model through SIR_model.R
##LOAD DATA
syn_data=readRDS("./syn_data.RDS")
y=syn_data$data$H
Rt=syn_data$data$Rt
hosp_prob=syn_data$hosp_prob

p=length(syn_data$theta_real)

runpCSMC_marginal=T
M = 10000 #Iterations
B = 300   #Particles per CSMC run

if(runpCSMC_marginal){
  resPCSMC_marginal =pSMCGibbs_marginal(y=y,R_t=Rt,hosp_prob=hosp_prob,p=p,B=B,M=M,
                      sseed = 3)
  saveRDS(resPCSMC_marginal,"SIR_PCSMC_marginal.rds")
}

#A couple of ideas: Implement RW version and check whether using only \thetastar_T and xstar_{t:T} is possible  
# (collapsed version where \theta_{1:(T-1)} are sampled at every iteration, also for the reference path)
