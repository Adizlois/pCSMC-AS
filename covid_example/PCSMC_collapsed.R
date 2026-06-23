# Author: Alfonso Diz-Lois Palomares
# Email: adpaloma@uio.no
library(extraDistr)
source("Aux_functions.R")
Rcpp::sourceCpp("get_cond.cpp")


#MAIN PARTIALLY COLLAPSED CSMC-AS ROUTINE

  col_pCSMC = function(y,xstar,thetastar,R_t,B,
                   hosp_prob,alpha=NULL,beta=NULL,AS=T)
  {
    nT = dim(xstar)[1]
    p = dim(xstar)[2]-1
    
    #Default priors hyperparameters
    if (is.null(alpha)) alpha <- 2
    if (is.null(beta))  beta  <- 4
    #Theta prior (gamma distribution)
    alphas=rep(alpha,p)
    betas=rep(beta,p)
    xsim = array(0,c(nT,p+1,B)) #Number of infections. Each "p" represents how many of the infections of that 
    #day come from each of the previous days. Last dimension is just the sum
    thetasim=array(NA,c(nT,p,B)) #Parameters sampled in the augmented space
    suffs=array(0,c(2,nT,p,B)) # 2 sufficients: counts and exposure (R_t*I_{t-m})
    #Set reference as first particle
    xsim[,,1] = xstar
    thetasim[,,1] =thetastar
    w = rep(0,B)
    #Seeding
    xsim[1,1:p,-1] = 0
    xsim[1,p+1,-1]=rpois(B-1,100)
    
    w = w+ dbinom(y[1],xsim[1,p+1,],hosp_prob,log=TRUE)
    w[is.na(w)]<-0
    ############
    #Thetas from prior for the first time point
    thetasim[1,,]<-rgamma(B*p,shape = alpha,rate = beta)
    #Be sure that if no samples are available during the first time points, 
    #the prior will be used
    for (j in 2:p){
      thetasim[j,,]<-thetasim[1,,]
    }
    thetasim[,p,]<-0.05 #Fixed theta to make the problem identifiable
    #Update sufficient
    suffs[1,1,,]<-alpha
    suffs[2,1,,]<-beta
    
    for(s in 2:nT)
    {
      #cat("\rTime:", s)
      #flush.console()
      a1 = 1
      #Ancestor sampling step
      if ((s>2)&(AS)){
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
        if (s<nT){
          kmax=min((p-1),s-1)
          wtilde=wtilde+get_cond_cpp(thetastar[1:kmax],alphas = suffs[1,s-1,1:kmax,],betas = suffs[2,(s-1),1:kmax,])
        }
        
        finite <- is.finite(wtilde)
        m <- max(wtilde[finite])              
        wtilde[finite] <- exp(wtilde[finite] - m)  
        wtilde[!finite]<-0
        a1 = sample(1:B,1,prob=wtilde)
      } ## END OF ANCESTOR SAMPLING STEP
      w = exp(w-max(w))
      
      #Resampling
      a = c(a1,sample(1:B,B-1,prob=w,replace=TRUE))
      xsim[1:(s-1),,] = xsim[1:(s-1),,a] 
      thetasim[1:(s-1),,]=thetasim[1:(s-1),,a]
      suffs[,1:(s-1),,]=suffs[,1:(s-1),,a]
      w = rep(0,B)
      #Forward model
      kmax <- min(p, s - 1)
      if (s>2){
        for(j in 1:kmax){
          if(j==p){
            thetasim[s,j,]=0.05 #Fixed theta
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
      #Update Sufficient statistic
      suffs[1,s,,]=suffs[1,s-1,,]+xsim[s,1:p,]
      suffs[2,s,1:kmax,]=suffs[2,s-1,1:kmax,]+R_t[s]*xsim[(s-1):(s-kmax),p+1,]
      if (kmax<p){
        suffs[2,s,(kmax+1):p,]=suffs[2,s-1,(kmax+1):p,]
      }  
      w = w+ dbinom(y[s],xsim[s,p+1,],hosp_prob,log=TRUE)
      w[is.infinite(w)]<--100 #Arbritrary
    }
    #Sample new reference at the last time point
    w = exp(w-max(w))
    k = sample(1:B,1,prob=w)
    return(list(latent=xsim[,,k],param=thetasim[nT,,k]))
  }

  
  #CSMC ROUTINE TO SAMPLE R_{1:T} GIVEN X_{1:T} AND THETA. WE ASSUME A RW
  #DYNAMIC MODEL IN R
  CSMC_R = function(x,Rstar,theta,B,sigma_R,
                    AS=T)
  {
    nT = length(Rstar)
    p = dim(x)[2]-1
    Rsim = array(NA,c(nT,B)) #Number of infections. Each "p" represents how many of the infections of that 
    #day come from each of the previous days. Last dimension is just the sum
    Rsim[,1] = Rstar
    w = rep(0,B)
    Rsim[1,-1] = exp(rnorm(B-1,0,sigma_R))
    for(s in 2:nT)
    {

      
      a1=1
      if (AS){ #ANCESTOR SAMPLING STEP
        wtilde = w
        wtilde=wtilde+dnorm(log(Rstar[s]),log(Rsim[s-1,]),sd=sigma_R,log=TRUE)
        wtilde[is.infinite(wtilde)]<--100 #Arbitrary
        wtilde <- exp(wtilde- max(wtilde))  
        a1 = sample(1:B,1,prob=wtilde)
      } #END OF ANCESTOR SAMPLING
      
      w = exp(w-max(w))
      a = c(a1,sample(1:B,B-1,prob=w,replace=TRUE))
      Rsim[1:(s-1),] = Rsim[1:(s-1),a] 
      w = rep(0,B)
      #Forward model
      Rsim[s,-1 ] <- exp(rnorm(B-1,log(Rsim[s-1,-1]),sd=sigma_R))
      
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
  thetasim = rgamma(p,2,40)
  thetasim[p]<-0.05
  thetaM[1,] = thetasim
  for(m in 1:M)
  {
    cat("\r Iteration ",m, " ")
    
    #Sample x and theta given R_{1:T} using col-pCSMC
    if(m==1){
    res<-col_pCSMC(y,xstar=xsim,thetastar=thetasim,R_t,B,
                     hosp_prob,alpha=alpha0,beta=beta0,AS=F)
    }
    else{
      res<-col_pCSMC(y,xstar=xsim,thetastar=thetasim,R_t,B,
                 hosp_prob,alpha=alpha0,beta=beta0,AS=T)
      
    }
    
    xsim=res$latent
    xsim[is.na(xsim)]<-0
    xsimM[m+1,,] = xsim
    thetasim=res$param
    thetaM[m+1,]=thetasim
    
    #Sample R_{1:T}|x and theta using CSMC
    if(m==1){
    R_t<-CSMC_R(x=xsim,Rstar=R_t,theta=thetasim,B=B,sigma_R=0.15,
                AS=F)
    
    }else{
      R_t<-CSMC_R(x=xsim,Rstar=R_t,theta=thetasim,B=B,sigma_R=0.15,
                  AS=T)
                }
    
    
    R_tM[m+1,]<-R_t
  }
  list(xsim=xsimM,thetaM=thetaM,R_tM=R_tM)
}

###### RUNNING THE MODEL..

library(data.table)
covid_data<-fread("incidence_national.csv")
covid_data<-covid_data[innleggelsesdato>as.Date("2021-02-01")&innleggelsesdato<as.Date("2021-03-15"),]
y         <- covid_data$N
hosp_prob <- 0.047
p         <- 8

#plot(y,ylab="Hospitalizations",xlab="Time")


##RUN PCSMC

  runPCSMC = T
  M = 90000 #Iterations
  B = 300  #Particles per iteration
  
  if(runPCSMC){
    resCSMC = SMCGibbs_R(y=y,hosp_prob=hosp_prob,p=p,B=B,M=M,
                         alpha0=2,beta0=40,seed_ = 132320)
    saveRDS(resCSMC,"SIR_PCSMC_covid.rds")
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