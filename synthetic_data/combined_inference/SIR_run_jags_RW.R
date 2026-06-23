library(rjags)
library(coda)
setwd("/home/gandalf/Documents/GitHubUiO/Branching_process/real_data/")
syn_data <- readRDS("./bastet_shorter/syn_data.RDS")
y         <- syn_data$data$H
hosp_prob <- syn_data$hosp_prob
p         <- length(syn_data$theta_real)

plot(y,ylab="Hospitalizations",xlab="Time")

data_list <- list(
  nT        = length(y),
  y         = as.integer(y),
  p         = p,
  hosp_prob = hosp_prob,
  sigma_R   = 0.25,
  log_R0     = 0,
  I0        = 2 * p,       # fixed seed
  alpha     = rep(2, p),   # for theta[2:p]
  beta      = rep(4, p)
)

inits_fun <- function() {
  # Initialize the free thetas (lags 2:p)
  th_free <- rgamma(p - 1, 1, 1)  # theta_free[1:(p-1)]
  
  p_clip  <- min(max(hosp_prob, 1e-6), 1 - 1e-6)
  I_guess <- (pmax(as.integer(round(y / p_clip)), y) * 2) + 1
  I_guess[1] <- NA  # I[1] is fixed to I0 in the model
  
  list(
    theta_free = th_free,
    I          = I_guess
    # Rt left for JAGS to initialize automatically
  )
}


params <- c("theta", "theta_free", "I", "Lambda", "Rt")

mod <- jags.model(
  file    = "SIR_model_jags_RW.jags",
  data    = data_list,
  inits   = inits_fun,
  n.chains = 3,
  n.adapt  = 20000
)

update(mod, 100000)  # burn-in

samps <- coda.samples(
  mod,
  variable.names = params,
  n.iter = 200000
  #thin   = 10
)

samps<-readRDS("SIR_samps_rw_long_proof.rds")
gelman.diag(samps,multivariate=FALSE)



# Posterior summaries
mcmc_list <- as.mcmc.list(samps)
nchains <- length(mcmc_list)
niter   <- nrow(mcmc_list[[1]])
npar    <- ncol(mcmc_list[[1]])

arr <- array(NA, dim = c(niter, npar, nchains))

for (i in 1:nchains) {
  arr[ , , i] <- as.matrix(mcmc_list[[i]])
}
saveRDS(samps,"./SIR_samps_rw_long_proof.rds")

par_names <- colnames(mcmc_list[[1]])

norm_R<-function(Rt,thetas){
  for (i in 1:dim(Rt)[1]){
    Rt[i,]<-Rt[i,]*sum(thetas[i,])
  }
  Rt
}
arr<-arr[,,1]
plot(apply(norm_R(arr[,300:449],arr[,450:453]),2,mean),ylab="Rt",xlab="t",type="l",col="blue")
lines(syn_data$data$Rt,col="black",lwd=1.5,lty=2)
legend("topright",legend=c("JAGS","real"),lty=c(1,2),lwd=c(1,1.2),col=c("blue","black"))



plot(apply(arr[,1:50],2,mean),type="l",col="blue")
lines(syn_data$data$new_infections,lty=2,lwd=1.2)

pcsmc<-readRDS("/home/gandalf/Downloads/SIR_PCSMC_fixed_both_csmc_1_25.rds")
lines(apply(norm_R(pcsmc$R_tM[15000:30000,],pcsmc$thetaM[15000:30000,]),2,mean),col="orange")
legend("topright",legend=c("pCSMC","JAGS","real"),lty=c(1,1,2),lwd=c(1,1,1.2),col=c("orange","blue","black"))


plot(density(arr[,451]/apply(arr[,450:453],1,sum)),main=bquote(theta[2]),col="blue")
abline(v=syn_data$theta_real[2],lwd=1.5,lty=2)
#lines(density(pcsmc$thetaM[15000:30000,4]/apply(pcsmc$thetaM[15000:30000,],1,sum)),col="orange")
legend("topright",legend=c("pCSMC","JAGS"),lty=c(1,1),lwd=c(1,1),col=c("orange","blue"))
