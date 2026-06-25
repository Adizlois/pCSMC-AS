library(extraDistr)
library(MASS)
library(tmvtnorm)
#Sample independent gamma and normalized
rindep_gamma <- function(alpha,beta) {
  y <- rgamma(length(alpha), shape = alpha, rate = beta)
  y / sum(y)
}
samptheta1<-function(xsim,R_t,alpha=NULL,beta=NULL){
  
  #Assuming independent gammas 
  p=dim(xsim)[2]-1
  nT=length(R_t)
  # Prior hyperparameters: α_m = τ/p by default, common rate β = 1e-6
  if (is.null(alpha)) alpha <- rep(2, p)
  if (is.null(beta))  beta  <- rep(4, p)
  
  newalpha=alpha+apply(xsim[2:nT,1:p],2,sum)
  #Exposures
  # Exposures by lag: E_m = sum_t R_t * X_{t-m}
  Xtot <- xsim[, p + 1]
  E <- numeric(p)
  for (m in 1:p) {
    E[m]=sum(R_t[(1+m):nT]*Xtot[1:(nT-m)])
  }
  newbeta=beta+E
  rgamma(length(alpha), shape = newalpha, rate = newbeta)
}
get_prior_sample <- function(incidence, prior) {
  n <- length(incidence)
  p <- length(prior)
  out_parts <- matrix(0, nrow = n, ncol = p)
  for (i in seq_len(n)) {
    sz <- incidence[i]
    if (sz == 0) {
      out_parts[i, ] <- rep(0, p)
    } else {
      out_parts[i, ] <- as.integer(rmultinom(1, size = sz, prob = prior)[, 1])
    }
  }
  
  out <- cbind(out_parts, total = as.integer(incidence))
  return(out)
}


get_cond<-function(thetastar,alphas,betas){
  B=dim(alphas)[2]
  mmat <- sapply(seq_len(B), function(b) {
    dgamma(thetastar,
           shape = alphas[, b],
           rate  = betas[, b], 
           log   = TRUE)
  })
  colSums(mmat)
}
  
#Function to sample proposals from the reference as in Finke&Thiery's i-RW-CSMC
sample_ref <- function(ref_sample, nsamples, l = 1) {
  p <- length(ref_sample)
  Z.0 = truncnorm::rtruncnorm(p,ref_sample,sqrt(l/(2*p)),a=0)
  output=array(NA,dim=c(nsamples,p))
  for(i in 1:p)
  {
    # then perturb each dimension independently  
    output[,i] = truncnorm::rtruncnorm(nsamples,Z.0[i],sqrt(l/(2*p)),a=0)
  }
  output
}

dmvnorm_vector_vs_matrix <- function(V, MU, Sigma, log = FALSE) {
  # V: p-vector
  # MU: p x B matrix of means
  # Sigma: p x p covariance matrix
  p <- length(V)
  B <- ncol(MU)
  
  # Cholesky decomposition 
  L <- chol(Sigma)
  log_det <- 2 * sum(log(diag(L)))
  L_inv <- backsolve(L, diag(p))
  Sigma_inv <- t(L_inv) %*% L_inv
  
  # Compute quadratic forms
  qf <- numeric(B)
  for (b in 1:B) {
    diff <- V - MU[,b]
    qf[b] <- sum((Sigma_inv %*% diff) * diff)
  }
  
  logdens <- -0.5 * (p*log(2*pi) + log_det + qf)
  if (log) return(logdens)
  return(exp(logdens))
}