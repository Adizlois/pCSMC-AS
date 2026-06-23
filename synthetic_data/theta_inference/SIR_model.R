simulate_SIR_with_kernel <- function(
    nT,                  # number of time steps (days)
    N,                  # total population size
    beta,               # transmission rate (per infectious "effective person" per day)
    infectious_period,  # mean infectious period (days); gamma = 1 / infectious_period
    hosp_prob,          # hospitalization probability
    theta,                  # infectivity kernel over 1:s days, nonnegative, sums to 1
    I0 = 10,            # initial infectious (seed as day-0 incidence)
    R0_init = 0,        # initial recovered
    stochastic = TRUE,  # use binomial stochastic transitions; if FALSE, deterministic
    seed = 1,           # RNG seed (used only if stochastic=TRUE)
    filename=NULL
) {
  stopifnot(nT >= 1, N > 0, beta >= 0, infectious_period > 0)
  theta <- as.numeric(theta)
  p <- length(theta)
  stopifnot(all(theta >= 0))
  #if (abs(sum(theta) - 1) > 1e-8) stop("w must sum to 1.")
  
  gamma <- 1 / infectious_period
  
  if (stochastic) set.seed(seed)
  
  # Storage
  X <- numeric(nT)     # new infections
  S <- numeric(nT)     # susceptibles
  I <- numeric(nT)     # infectious
  R <- numeric(nT)     # removed
  J <- numeric(nT)     # kernel-weighted infectiousness (sum w_m * X_{t-m})
  H <- numeric(nT)
  Rt<- numeric(nT)    # beta * infectious_period * S_{t-1}/N
  
  # Initialize
  # Seed day-0 incidence to bootstrap the kernel
  past_inc <- rep(0, p)
  past_inc[1] <- I0
  
  S_prev <- N - I0 - R0_init
  I_prev <- I0
  R_prev <- R0_init
  
  for (t in seq_len(nT)) {
    # Kernel-weighted infectiousness from past incidence
    J[t] <- sum(theta * past_inc)
    
    # Effective reproduction at time t (SIR)
    Rt[t] <- beta[t] * S_prev / N
    
    # Incidence update
    if (stochastic) {
      p_inf <- 1 - exp(-beta[t] * J[t] / N)
      p_inf <- max(0, min(1, p_inf))
      X[t] <- rbinom(1, size = round(S_prev), prob = p_inf)
    } else {
      X[t] <- beta[t] * (S_prev / N) * J[t]
    }
    
    # Recoveries
    if (stochastic) {
      p_rec <- 1 - exp(-gamma)  # per-day removal probability
      p_rec <- max(0, min(1, p_rec))
      Y_t <- rbinom(1, size = round(I_prev), prob = p_rec)
    } else {
      Y_t <- gamma * I_prev
    }
    
    # Compartment updates
    S[t] <- max(0, S_prev - X[t])
    I[t] <- max(0, I_prev + X[t] - Y_t)
    H[t] <- rbinom(1,X[t],prob = hosp_prob)
    R[t] <- min(N, R_prev + Y_t)
    
    
    # Roll past incidence buffer with today's X[t]
    if (p > 1) {
      past_inc <- c(X[t], head(past_inc, p - 1))
    } else {
      past_inc[1] <- X[t]
    }
    
    # Prepare for next day
    S_prev <- S[t]
    I_prev <- I[t]
    R_prev <- R[t]
  }
  
  res<-list(data=data.frame(
    t = 1:nT,
    new_infections = X,
    J = J,
    H = H,
    S = S,
    I = I,
    R = R,
    Rt = Rt),
    theta_real=theta,
    hosp_prob=hosp_prob,
    infectious_period=infectious_period,
    I0=I0
  )
  if (!is.null(filename))
    saveRDS(object = res,file = filename)
  return(res)
}

###GENERATE DATA 

N <- 100000 # Population
nT <- 50; sigma <- 0.1 ; phi=0.95 # Using an AR1 model to generate the betas
set.seed(411)  #SEEDING

log_beta <- as.numeric(arima.sim(
  model = list(ar = phi),
  n = nT,
  rand.gen = function(n) rnorm(n, mean = 0, sd = sigma)
)) 
#MODEL TO GENERATE THE BETAS IN LOG SCALE
beta <- exp(log_beta) #TO EXP SCALE
beta=beta+0.9 #JUST A BIT OF INCREASE TO THE SIGNAL so that the outbreak has a reasonable size

infectious_period <- 3 # (not really important in this case)
p <- 2 #NUMBER OF DIMENSIONS CONSIDERED


#Real parameters
theta_real <- rgamma(p,2,4) 
#We normalize
theta_real=theta_real/sum(theta_real)

hosp_prob <- 0.1 #Hospitalization probability
sim <- simulate_SIR_with_kernel(
   nT = length(beta), N = N, beta = beta, infectious_period = infectious_period,hosp_prob=hosp_prob,
   theta = theta_real, I0 = 10, R0_init = 0, stochastic = TRUE, seed = 222,
  filename="syn_data.RDS"
 )
 plot(sim$data$H)


