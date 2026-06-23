setwd("/mn/sarpanitu/ansatte-u2/adpaloma/Documents/branching_real/")
source("SIR_model.R")
library(extraDistr)
source("Aux_functions.R")
Rcpp::sourceCpp("get_cond.cpp")


#MAIN CSMC ROUTINE

  pCSMC = function(y,xstar,thetastar,R_t,B,
                   hosp_prob,alpha=NULL,beta=NULL,AS=T)
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
    thetasim=array(NA,c(nT,p,B))
    suffs=array(0,c(2,nT,p,B)) # 2 sufficients: counts and exposure (R_t*I_{t-m})
    #day come from each of the previous days. Last dimension is just the sum
    xsim[,,1] = xstar
    thetasim[,,1] =thetastar
    w = rep(0,B)
    #Seeding
    xsim[1,1:p,-1] = 2
    xsim[1,p+1,-1]=2*p
    
    w = w+ dbinom(y[1],xsim[1,p+1,],hosp_prob,log=TRUE)
    w[is.na(w)]<-0
    #For the first p points, we will just use the prior two options: based on the reference (as in Finke%Thiery) or pure prior
    ############
    #Based on the reference
    #thetasim[1,,-1]<-sample_ref(thetastar[nT,],nsamples=B-1,l=0.1)
    
    ############
    thetasim[1,,]<-rgamma(B*p,shape = alpha,rate = beta)
    
    for (j in 2:p){
      thetasim[j,,]<-thetasim[1,,]
    }
    thetasim[,1,]<-1
    #Update sufficient
    suffs[1,1,,]<-alpha
    suffs[2,1,,]<-beta
    
    for(s in 2:nT)
    {
      #cat("\rTime:", s)
      flush.console()
      a1 = 1
      #Ancestor sampling step
      #if ((s>(p+1))){#&(s%%2==0)){
      if ((s>2)&(AS)){#&(s%%2==0)){
        wtilde = w
        #For p(xstar[t,]|xsim[1:t-1,,],\thetasta[nT,])
        for (m in 1:min(s-1,p)){
          wtilde=wtilde+dpois(xstar[s,m],R_t[s]*thetastar[m]*xsim[s-m,p+1,],log=TRUE)
        }
        
        #The other terms p(xstar[t+1:t+p-1,]|) We only need the terms that depend on the particle!
        for (j in 2:p){
          if ((s+j-1)<=nT){
            kmax <- min(p-j+1, s - 1)
            for (m in 1:kmax){
              wtilde=wtilde+dpois(xstar[s+j-1,j-1+m],R_t[s+j-1]*thetastar[j-1+m]*xsim[s-m,p+1,],log=TRUE)
            }
          }
        }
        
        
        #Theta terms
        #First the one that corresponds to p(\theta_T|x_{1:t-1}^i)
        if (s<nT){
          kmax=min(p,s-1)
          #if(kmax>2){
          wtilde=wtilde+get_cond_cpp(thetastar[1:kmax],alphas = suffs[1,s-1,1:kmax,],betas = suffs[2,(s-1),1:kmax,])
          #}
        }
        
        finite <- is.finite(wtilde)
        m <- max(wtilde[finite])              
        wtilde[finite] <- exp(wtilde[finite] - m)  
        wtilde[!finite]<-0
        a1 = sample(1:B,1,prob=wtilde)
      }
      
      #Now we also need the terms that correspond to the probability of the thetas in the reference
      #given the particle for t:nT
      w = exp(w-max(w))
      a = c(a1,sample(1:B,B-1,prob=w,replace=TRUE))
      xsim[1:(s-1),,] = xsim[1:(s-1),,a] 
      thetasim[1:(s-1),,]=thetasim[1:(s-1),,a]
      suffs[,1:(s-1),,]=suffs[,1:(s-1),,a]
      w = rep(0,B)
      #Forward model
      
      kmax <- min(p, s - 1)
      
      
      if (s>2){
        for(j in 1:kmax){
          if(j==1){
            thetasim[s,j,]=1
          }else{
          thetasim[s,j,]=rgamma(B,shape = suffs[1,s-1,j,],rate=suffs[2,s-1,j,])
          }
        }
      }
      for(j in 1:kmax){
        xsim[s,j,-1 ] <- rpois(B-1,R_t[s]*thetasim[s,j,-1]*xsim[s-j,p+1,-1])
      }
      xsim[s:nT,,1]<-xstar[s:nT,]
      thetasim[nT,,1]=thetastar
      xsim[s, ,-1][is.na(xsim[s, ,-1])] <- 0
      xsim[s,p+1,-1]=apply(xsim[s,1:p,-1],2,sum,na.rm=T)
      suffs[1,s,,]=suffs[1,s-1,,]+xsim[s,1:p,]
      suffs[2,s,1:kmax,]=suffs[2,s-1,1:kmax,]+R_t[s]*xsim[(s-1):(s-kmax),p+1,]
      if (kmax<p){
        suffs[2,s,(kmax+1):p,]=suffs[2,s-1,(kmax+1):p,]
      }  
      
      w = w+ dbinom(y[s],xsim[s,p+1,],hosp_prob,log=TRUE)
      w[is.infinite(w)]<--100 #Arbritrary
    }
    w = exp(w-max(w))
    k = sample(1:B,1,prob=w)
    return(list(latent=xsim[,,k],param=thetasim[nT,,k]))
  }

  
  #MAIN CSMC ROUTINE
  CSMC_R = function(x,Rstar,theta,B,sigma_R,
                    AS=T)
  {
    nT = length(Rstar)
    p = dim(x)[2]-1
    Rsim = array(NA,c(nT,B)) #Number of infections. Each "p" represents how many of the infections of that 
    #day come from each of the previous days. Last dimension is just the sum
    Rsim[,1] = Rstar
    w = rep(0,B)
    Rsim[1,-1] = runif(B-1,0.5,2)
    for(s in 2:nT)
    {
      #cat("\rTime:", s)
      #flush.console()
      a1=1
      if (AS){
        wtilde = w
        kmax <- min(p, s - 1)
        #For p(xstar[t,]|xsim[1:t-1,,])
        wtilde=wtilde+dnorm(log(Rstar[s]),0.95*log(Rsim[s-1,]),sd=sigma_R,log=TRUE)
        
        wtilde[is.infinite(wtilde)]<--100 #Arbitrary
        
        wtilde <- exp(wtilde- max(wtilde))  
        a1 = sample(1:B,1,prob=wtilde)
      }
      
      w = exp(w-max(w))
      a = c(a1,sample(1:B,B-1,prob=w,replace=TRUE))
      Rsim[1:(s-1),] = Rsim[1:(s-1),a] 
      w = rep(0,B)
      #Forward model
      Rsim[s,-1 ] <- exp(rnorm(B-1,0.95*log(Rsim[s-1,-1]),sd=sigma_R))
      
      kmax <- min(p, s - 1)
      
      for(k in 1:kmax){
        w = w+ dpois(x[s,k],Rsim[s,]*theta[k]*x[(s-k),p+1],log=TRUE)
        w[is.infinite(w)]<--100 #Arbitrary
      }
    }
    w = exp(w-max(w))
    k = sample(1:B,1,prob=w)
    Rsim[,k]
  }
  
  
  
##MAIN FUNCTION


SMCGibbs_R = function(y,hosp_prob,p,alpha0=2,beta0=4,B=100,M=100,seed_=100)
{
  set.seed(seed_)
  nT = length(y)
  xsimM = array(NA,c(M+1,nT,p+1))
  thetaM = matrix(NA,nrow=M+1,ncol=p)
  R_tM=matrix(NA,nrow=M+1,ncol=nT)
  # Initializing around observation
  xsim = get_prior_sample(round(y/hosp_prob),prior=rep(1,p)) 
  R_t=c(0.5,(y[2:nT]-y[1:(nT-1)])/y[1:(nT-1)])
  R_t[(is.na(R_t))|R_t<=0|is.infinite(R_t)]<-0.5
  
  R_tM[1,]<-R_t
  xsimM[1,,]<-xsim 
  thetasim = rgamma(p,2,4)
  thetasim[1]<-1
  thetaM[1,] = thetasim
  for(m in 1:M)
  {
    cat("\r Iteration ",m, " ")
    
    #Normalization
    #thetasim=thetasim/sum(thetasim)
    
    #Sample x using CSMC
    if(m==1){
    res<-pCSMC(y,xstar=xsim,thetastar=thetasim,R_t,B,
                     hosp_prob,alpha=NULL,beta=NULL,AS=F)
    }
    else{
      res<-pCSMC(y,xstar=xsim,thetastar=thetasim,R_t,B,
                 hosp_prob,alpha=NULL,beta=NULL,AS=T)
      
    }
    
    xsim=res$latent
    xsim[is.na(xsim)]<-0
    xsimM[m+1,,] = xsim
    thetasim=res$param
    thetaM[m+1,]=thetasim
    if(m==1){
    R_t<-CSMC_R(x=xsim,Rstar=R_t,theta=thetasim,B=B,sigma_R=0.25,
                AS=F)
    
    }else{
      R_t<-CSMC_R(x=xsim,Rstar=R_t,theta=thetasim,B=B,sigma_R=0.25,
                  AS=T)
                }
    
    
    #Sample theta|xsim,R_t
    
    #Some plots along the way
    # if (m%%100==0){
    #   flush.console()
    #   plot_param_traces(thetaM[1:m,]/apply(thetaM[1:m,],1,sum),sim$theta_real,trueI=sim$data$new_infections,simI=xsim[,p+1])
    # }
    R_tM[m+1,]<-R_t
  }
  list(xsim=xsimM,thetaM=thetaM,R_tM=R_tM)
}

# Example usage


###GENERATE DATA (ALTERNATIVELY ONE COULD READ THE "syn_data.RDS" FILE DIRECTLY)

N <- 1000000 # Population
nT <- 70; sigma <- 0.02 ; phi=0.95 # NUMBER OF DAYS AND DATA TO GENERATE THE BETAS (RT) WITH AN AUTOREGRESSIVE APPROACH
#set.seed(911)  #SEEDING
set.seed(411)  #SEEDING

log_beta<-runif(1,0.5,0.8)
for (j in 2:nT){
  log_beta<-c(log_beta,sample(c(phi,1.1),size = 1,prob = c(0.8,0.2))*log_beta[length(log_beta)]+rnorm(1,sd=sigma))
}
plot(exp(log_beta))

#log_beta <- as.numeric(arima.sim(
#  model = list(ar = phi),
#  n = nT,
#  rand.gen = function(n) rnorm(n, mean = 0, sd = sigma)
#)) #MODEL TO GENERATE THE BETAS IN LOG SCALE
beta <- exp(log_beta) #TO EXP SCALE
#beta=beta+0.9 #JUST A BIT OF INCREASE TO THE SIGNAL

infectious_period <- 3 # (not really important)
p <- 4 #NUMBER OF DIMENSIONS CONSIDERED

theta_real <- rgamma(p,2,4) #Real parameters
theta_real=theta_real/sum(theta_real)

hosp_prob <- 0.1 #Hospitalization probability
sim <- simulate_SIR_with_kernel(
  nT = length(beta), N = N, beta = beta, infectious_period = infectious_period,hosp_prob=hosp_prob,
  theta = theta_real, I0 = 10, R0_init = 0, stochastic = TRUE, seed = 222,#323
  filename="syn_data.RDS"
)
plot(sim$data$H)
plot(sim$data$Rt)
sim=readRDS("./syn_data.RDS")
y=sim$data$H
Rt=sim$data$Rt
hosp_prob=sim$hosp_prob
theta_real=sim$theta_real

# p=length(sim$theta_real)

##RUN PGAS

  runPCSMC = T
  M = 60000
  B = 300
  
  if(runPCSMC){
    resCSMC = SMCGibbs_R(y=sim$data$H,hosp_prob=sim$hosp_prob,p=length(sim$theta_real),B=B,M=M)
    saveRDS(resCSMC,"SIR_PCSMC_fixed_both_csmc_1_long.rds")
  }
  
  
  # 
  # burnin=1000
  # # summary statistics over M for each t
  # 
  # R_mean <- apply(resPGAS$R_tM[burnin:M,]/rowSums(resPGAS$thetaM[burnin:M,]), 2, mean)
  # R_low  <- apply(resPGAS$R_tM[burnin:M,]/rowSums(resPGAS$thetaM[burnin:M,]), 2, quantile, probs = 0.025)
  # R_high <- apply(resPGAS$R_tM[burnin:M,]/rowSums(resPGAS$thetaM[burnin:M,]), 2, quantile, probs = 0.975)
  # 
  # t_vec <- 1:ncol(resPGAS$R_tM[burnin:M,])
  # 
  # # base plot
  # plot(t_vec, R_mean, type = "l", lwd = 2,
  #      xlab = "t", ylab = "R_t",
  #      ylim = range(c(R_low, R_high)))
  # 
  # # add 95% intervals as lines
  # lines(t_vec, R_low,  lty = 2)
  # lines(t_vec, R_high, lty = 2)
  # 
  # lines(sim$data$Rt,lty=2,col="blue",lwd=1.1)
  # 
  # burnin=500
  # # summary statistics over M for each t
  # R_mean <- apply(resPGAS$R_tM[burnin:M,], 2, mean)
  # R_low  <- apply(resPGAS$R_tM[burnin:M,], 2, quantile, probs = 0.025)
  # R_high <- apply(resPGAS$R_tM[burnin:M,], 2, quantile, probs = 0.975)
  # 
  # t_vec <- 1:ncol(resPGAS$R_tM[burnin:M,])
  # 
  # # base plot
  # plot(t_vec, R_mean, type = "l", lwd = 2,
  #      xlab = "t", ylab = "R_t",
  #      ylim = range(c(R_low, R_high)))
  # 
  # # add 95% intervals as lines
  # lines(t_vec, R_low,  lty = 2)
  # lines(t_vec, R_high, lty = 2)
  # 
  # lines(sim$data$Rt,lty=2,col="blue",lwd=1.1)
  # 
  # #Thetas
  # par(mfrow=c(round(p/2),2))
  # for (j in 1:p){
  #   plot(density(resPGAS$thetaM[burnin:M,j]/apply(resPGAS$thetaM[burnin:M,],1,sum)),col="blue",
  #        main= bquote(theta[.(j)]))
  # #  plot(density(resPGAS$thetaM[burnin:M,j]),col="blue",
  # #              main= bquote(theta[.(j)]))
  #   abline(v=theta_real[j],lty=2,lwd=1.2)
  # }
  # layout(1)
  # 
