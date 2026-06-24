# Author: Alfonso Diz-Lois Palomares
# Email: adpaloma@uio.no
library(rjags)
library(coda)

#The script assumes the synthetic data is available as syn_data.RDS. 
#Otherwise, run the model through SIR_model.R
syn_data <- readRDS("syn_data.RDS")
y         <- syn_data$data$H
hosp_prob <- syn_data$hosp_prob
p         <- length(syn_data$theta_real)

# Inputs 
# y: integer vector length T with hospitalizations
# p: kernel length (number of lags)
# hosp_prob: scalar in [0, 1]
# sigma_R : SD of RW for R_t in the log scale
# log_R0 : mean at time 1
# I0 : seeding
# alpha : hyperparameter prior for theta (Gamma distributed)
# beta : hyperparameter prior for theta (Gamma distributed)

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

#Initialize
mod <- jags.model(
  file    = "SIR_model_jags_RW.jags",
  data    = data_list,
  inits   = inits_fun,
  n.chains = 3,
  n.adapt  = 20000
)
#Burn-in
update(mod, 100000)  

#Sample
samps <- coda.samples(
  mod,
  variable.names = params,
  n.iter = 200000
)

saveRDS(samps,"./Jags.rds")
