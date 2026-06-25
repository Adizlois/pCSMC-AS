# Author: Alfonso Diz-Lois Palomares
# Email: adpaloma@uio.no

source("Aux_functions.R")
Rcpp::sourceCpp("get_cond.cpp")
library(matrixStats)
###IMPLEMENTATION OF PCSMC WITH ANCESTOR SAMPLING
pCSMC = function(y,xstar,thetastar,R_t,B,
                 hosp_prob,alpha=NULL,beta=NULL)
{
  nT = dim(xstar)[1]
  p = dim(xstar)[2]-1
  #Default priors
  if (is.null(alpha)) alpha <- 2
  if (is.null(beta))  beta  <- 4
  #Theta prior (gamma distribution)
  alphas=rep(alpha,p)
  betas=rep(beta,p)
  #Objects
  xsim = array(0,c(nT,p+1,B)) #Number of infections. Each "p" represents how many of the infections of that 
  #day come from each of the previous days. Last dimension (p+1) is just the sum 
  thetasim=array(NA,c(nT,p,B)) #Parameters (augmented space)
  suffs=array(0,c(2,nT,p,B)) # 2 sufficients: counts and exposure (R_t*I_{t-m})
  #Set the reference path as particle 1
  xsim[,,1] = xstar
  thetasim[,,1] =thetastar
  #Weights
  w = rep(0,B)
  #Seeding
  xsim[1,1:p,-1] = 2
  xsim[1,p+1,-1]=2*p
  
  w = w+ dbinom(y[1],xsim[1,p+1,],hosp_prob,log=TRUE)
  w[is.na(w)]<-0

  ############
  thetasim[1,,-1]<-rgamma((B-1)*p,shape = alpha,rate = beta) #Initial values
  
  #Be sure that if no samples are available during the first time points, 
  #the prior will be used
  for (j in 2:p){
    thetasim[j,,-1]<-thetasim[1,,-1]
  }
  #Update sufficient
  suffs[1,1,,]<-alpha
  suffs[2,1,,]<-beta
  
  for(s in 2:nT)
  {
    #cat("\rTime:", s)
    #flush.console()
    a1 = 1
    
    #Ancestor sampling step
    if (s>2){
    wtilde = w
    #For p(xstar[t,]|xsim[1:t-1,,],\thetasta[nT,])
    for (m in 1:min(s-1,p)){
      wtilde=wtilde+dpois(xstar[s,m],R_t[s]*thetastar[nT,m]*xsim[s-m,p+1,],log=TRUE)
    }

    #The other terms p(xstar[t+1:t+p-1,]|) We only need the terms that depend on the particle!
    for (j in 2:p){
      if ((s+j-1)<=nT){
        kmax <- min(p-j+1, s - 1)
        for (m in 1:kmax){
          wtilde=wtilde+dpois(xstar[s+j-1,j-1+m],R_t[s+j-1]*thetastar[nT,j-1+m]*xsim[s-m,p+1,],log=TRUE)
        }
      }
    }

    #Theta terms
    #First the one that corresponds to p(\theta_T|x_{1:t-1}^i)
    if (s<nT){
        kmax=min(p,s-1)
       wtilde=wtilde+get_cond_cpp(thetastar[nT,1:kmax],alphas = suffs[1,s-1,1:kmax,],betas = suffs[2,s-1,1:kmax,])

    #Now the other estimates
    #First the cumulative sum of xstar to get the alphas
       combined_particle<-xsim
       combined_particle[s:nT, , ] <- xstar[s:nT,]
       allalphas=alpha+apply(combined_particle[2:nT, 1:p, , drop = FALSE],
                       MARGIN = c(2, 3),
                       FUN    = cumsum)

       for (ss in s:(nT-1)){
        alphas=allalphas[ss-1,,]
        betas <- matrix(0, p, B)
        for (j in 1:p) {
          betas[j, ] <- beta+crossprod(R_t[(j+1):ss],xsim[1:(ss-j),p+1, ])
        }
        wtilde=wtilde+get_cond_cpp(thetastar[ss,],alphas=alphas,betas=betas)
      }
    }

    finite <- is.finite(wtilde)
    m <- max(wtilde[finite])              
    wtilde[finite] <- exp(wtilde[finite] - m)  
    wtilde[!finite]<-0
    a1 = sample(1:B,1,prob=wtilde)
    } #End of Ancestor Sampling step

    
    w = exp(w-max(w))
    ##Resampling step
    a = c(a1,sample(1:B,B-1,prob=w,replace=TRUE))
    xsim[1:(s-1),,] = xsim[1:(s-1),,a] 
    thetasim[1:(s-1),,]=thetasim[1:(s-1),,a]
    suffs[,1:(s-1),,]=suffs[,1:(s-1),,a]
    ##
    w = rep(0,B)
    #Forward model
    kmax <- min(p, s - 1)
    if (s>2){
      for(j in 1:kmax){
        thetasim[s,j,-1]=rgamma(B-1,shape = suffs[1,s-1,j,-1],rate=suffs[2,s-1,j,-1])
      }
    }
    for(j in 1:kmax){
      xsim[s,j,-1 ] <- rpois(B-1,R_t[s]*thetasim[s,j,-1]*xsim[s-j,p+1,-1])
    }
    xsim[s:nT,,1]<-xstar[s:nT,]
    
    xsim[s, ,-1][is.na(xsim[s, ,-1])] <- 0
    xsim[s,p+1,-1]=apply(xsim[s,1:p,-1],2,sum,na.rm=T)
    
    #Update Sufficient
    suffs[1,s,,]=suffs[1,s-1,,]+xsim[s,1:p,]
    suffs[2,s,1:kmax,]=suffs[2,s-1,1:kmax,]+R_t[s]*xsim[(s-1):(s-kmax),p+1,]
    if (kmax<p){
      suffs[2,s,(kmax+1):p,]=suffs[2,s-1,(kmax+1):p,]
    }  
    
    w = w+ dbinom(y[s],xsim[s,p+1,],hosp_prob,log=TRUE)
    w[is.infinite(w)]<--100 #Arbritrary
  }
  #Resampling on the last time step
  w = exp(w-max(w))
  k = sample(1:B,1,prob=w)
  return(list(latent=xsim[,,k],param=thetasim[,,k]))
}

#MAIN FUNCTION
pSMCGibbs = function(y,R_t,hosp_prob,p,B=100,M=100,sseed=10)
{
  set.seed(sseed)
  nT = length(y)
  xsimM = array(NA,c(M+1,nT,p+1))
  thetaM = array(NA,c(M+1,nT,p))
  # Initializing xsim around observation
  xsim = get_prior_sample(round(y/hosp_prob),prior=rep(1,p)) 
  xsim[1,1:p] = 2
  xsim[1,p+1]=2*p
  xsimM[1,,]<-xsim 
  
  #Initial random theta
  thetaM[1,,]<-rgamma(p,1,1)
  
  #Iterations loop
  for(m in 1:M)
  {
    cat("\r",m, " ")
    #Sample x using CSMC
    res = pCSMC(y,xstar=xsim,thetastar=thetaM[m,,],R_t,B,
                hosp_prob)
    xsim=res$latent
    xsim[is.na(xsim)]<-0
    xsimM[m+1,,] = xsim
    thetaM[m+1,,]=res$param
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

runpCSMC=T

M = 10000 #Iterations
B = 300  #Samples per iteration

if(runpCSMC){
  resPCSMC =pSMCGibbs(y=y,R_t=Rt,hosp_prob=hosp_prob,p=p,B=B,M=M,
                      sseed = 3)
  saveRDS(resPCSMC,"SIR_PCSMC.rds")
}
