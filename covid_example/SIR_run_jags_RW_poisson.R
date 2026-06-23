# Author: Alfonso Diz-Lois Palomares
# Email: adpaloma@uio.no
library(rjags)
library(coda)
library(data.table)

#Loading real data
covid_data<-fread("incidence_national.csv")
covid_data<-covid_data[innleggelsesdato>as.Date("2021-02-01")&innleggelsesdato<as.Date("2021-03-15"),]
y         <- covid_data$N

hosp_prob <- 0.047 #Hospitalization probability
p         <- 8 #Number of lags considered

#Input data
data_list <- list(
  nT        = length(y),
  y         = as.integer(y),
  p         = p,
  hosp_prob = hosp_prob,
  sigma_R   = 0.15, #Random walk SD for log R_t
  log_R0     = 0, #Initial value mean (R1=RW(log_R0,sigma_R))
  alpha     = rep(2, p),   # prior hyperparameters
  beta      = rep(40,p)
)

set.seed(102)

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

#Latent & Parameters
params <- c("theta", "theta_free", "I", "Lambda", "Rt")

mod <- jags.model(
  file    = "SIR_model_jags_RW_poisson.jags",
  data    = data_list,
  inits   = inits_fun,
  n.chains = 3,
  n.adapt  = 20000
)

update(mod, 60000)  # burn-in

samps <- coda.samples(
  mod,
  variable.names = params,
  n.iter = 70000,
  thin=3
)
saveRDS(samps, file = "jags_samples.rds")
