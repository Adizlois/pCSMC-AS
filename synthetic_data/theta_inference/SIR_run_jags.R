# Author: Alfonso Diz-Lois Palomares
# Email: adpaloma@uio.no
library(rjags)
library(coda)

# Inputs 
# y: integer vector length T with hospitalizations
# Rt: numeric vector length T with R_t
# p: kernel length (number of lags)
# hosp_prob: scalar in [0, 1]

#The script assumes the synthetic data is available as syn_data.RDS. 
#Otherwise, run the model through SIR_model.R
syn_data=readRDS("./syn_data.RDS")
y=syn_data$data$H
Rt=syn_data$data$Rt
hosp_prob=syn_data$hosp_prob
p=length(syn_data$theta_real)

data_list <- list(
  nT =length(y),
  y = as.integer(y),
  Rt = as.numeric(Rt),
  p = p,
  hosp_prob = hosp_prob,
  I0 = 2*p,                 # fixed seed
  alpha=rep(2,p),
  beta=rep(4,p) 
)



inits_fun <- function() {
  th <- rgamma(p,1,1)
  p_clip <- min(max(hosp_prob, 1e-6), 1 - 1e-6)
  I_guess <- (pmax(as.integer(round(y / p_clip)), y)*2)+1
  I_guess[1]<-NA
  init <- list(theta = th,
               I=I_guess)
  init
}

params <- c("theta", "I", "Lambda")

mod <- jags.model(
  file = "SIR_model_jags.jags",
  data = data_list,
  inits = inits_fun,
  n.chains = 3,
  n.adapt = 2000
)

update(mod, 50000)  # burn-in

samps <- coda.samples(
  mod,
  variable.names = params,
  n.iter = 50000,
  thin = 3
)
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
saveRDS(arr,"./SIR_jags.rds")
