# Author: Alfonso Diz-Lois Palomares
# Email: adpaloma@uio.no
library(extraDistr)
source("Aux_functions.R")


#CSMC ALGORITHM
CSMC = function(y,xstar,theta,R_t,B,
                hosp_prob,AS=T)
{
  nT = dim(xstar)[1]
  p = dim(xstar)[2]-1
  xsim = array(NA,c(nT,p+1,B)) #Number of infections. Each "p" represents how many of the infections of that 
  #day come from each ofthe previous days. Last dimension is just the sum
  xsim[,,1] = xstar #Set reference as particle 1
  w = rep(0,B)
  #Seeding
  xsim[1,1:p,-1] = 2
  xsim[1,p+1,-1]=2*p
  w = w+ dbinom(y[1],xsim[1,p+1,],hosp_prob,log=TRUE)
  for(s in 2:nT)
  {
    #cat("\rTime:", s)
    #flush.console()
    a1=1
    
    #Ancestor sampling step
    if (AS){
    wtilde = w
    kmax <- min(p, s - 1)
    #For p(xstar[t,]|xsim[1:t-1,,])
    for (m in 1:kmax){
      wtilde=wtilde+dpois(xstar[s,m],R_t[s]*theta[m]*xsim[s-m,p+1,],log=TRUE)
    }
    
    #The other terms p(xstar[t+1:t+p-1,]|) We only need the terms that depend on the particle!
    for (j in 2:p){
      if ((s+j-1)<=nT){
        kmax <- min(p-j+1, s - 1)
        for (m in 1:kmax){
          wtilde=wtilde+dpois(xstar[s+j-1,j-1+m],R_t[s+j-1]*theta[j-1+m]*xsim[s-m,p+1,],log=TRUE)
        }
      }
    }
    finite <- is.finite(wtilde)
    m <- max(wtilde[finite])              
    wtilde[finite] <- exp(wtilde[finite] - m)  
    wtilde[!finite]<-0
    a1 = sample(1:B,1,prob=wtilde)
    } ##End of Ancestor sampling step

    w = exp(w-max(w))
    #Resampling
    a = c(a1,sample(1:B,B-1,prob=w,replace=TRUE))
    xsim[1:(s-1),,] = xsim[1:(s-1),,a] 
    w = rep(0,B)
    #Forward model
    kmax <- min(p, s - 1)
    for(j in 1:kmax){
      xsim[s,j,-1 ] <- rpois(B-1,R_t[s]*theta[j]*xsim[s-j,p+1,-1])
    }
    xsim[s,p+1,-1]=apply(xsim[s,1:p,-1],2,sum,na.rm=T)
    w = w+ dbinom(y[s],xsim[s,p+1,],hosp_prob,log=TRUE)
  }
  #Sampling the new reference at the last time point
  w = exp(w-max(w))
  k = sample(1:B,1,prob=w)
  xsim[,,k]
}

##MAIN FUNCTION
SMCGibbs = function(y,R_t,hosp_prob,p,B=100,M=100,seed_=121)
{
  set.seed(seed_)
  nT = length(y)
  xsimM = array(NA,c(M,nT,p+1))
  thetaM = matrix(NA,nrow=M,ncol=p)
  # Initializing around observation
  xsim = get_prior_sample(round(y/hosp_prob),prior=rep(1,p)) 
  xsimM[1,,]<-xsim 
  #First theta sample
  thetasim = rgamma(p,2,4)
  
  for(m in 1:M)
  {
    cat("\r Iteration ",m, " ")
    #Sample theta from conditional posterior
    thetaM[m,] = thetasim
    
    
    #Sample x using CSMC (first iteration without ancestor sampling to avoid
    #potential incompatibilities between the first reference and the samples)
    if (m==1){
    xsim = CSMC(y,xstar=xsim,theta=thetasim,R_t,B,
                hosp_prob,AS=F)
    }else{
      xsim = CSMC(y,xstar=xsim,theta=thetasim,R_t,B,
                  hosp_prob,AS=T)
    }
    
    xsimM[m,,] = xsim
    #Sample theta from the conditional distribution (function from aux_functions.R)
    thetasim=samptheta1(xsim = xsim,R_t=R_t)
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

##RUN PGAS

  runPGAS = T
  M = 40000 #Iterations
  B = 300 #Particles
  
  if(runPGAS){
    resPGAS = SMCGibbs(y=y,R_t=Rt,hosp_prob=hosp_prob,p=p,B=B,M=M,seed_=132)
    saveRDS(resPGAS,"SIR_PGAS.rds")
  }
