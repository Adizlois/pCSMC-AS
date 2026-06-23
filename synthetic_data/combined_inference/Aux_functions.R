library(extraDistr)
library(MASS)
library(tmvtnorm)
source("SIR_model.R")

plot_param_traces <- function(theta_mat, theta_true, simI,trueI,param_names = NULL,
                              main = "Parameter trace plots",
                              col_trace = "#1f77b4", col_true = "#d62728",
                              lwd_trace = 1.5, lty_true = 2,add_I=T) {
  # Input checks
  if (!is.matrix(theta_mat)) stop("theta_mat must be a matrix (M x p).")
  M <- nrow(theta_mat)
  p <- ncol(theta_mat)
  
  if (length(theta_true) != p) stop("theta_true must have length equal to ncol(theta_mat).")
  
  # Names
  if (is.null(param_names)) {
    param_names <- colnames(theta_mat)
    if (is.null(param_names)) {
      param_names <- paste0("param", seq_len(p))
    }
  } else if (length(param_names) != p) {
    stop("param_names must be NULL or length equal to ncol(theta_mat).")
  }
  
  # Layout: try to make it roughly square
  ncol_plot <- 2#ceiling(sqrt(p+1))
  nrow_plot <- ceiling((p+1)/2)#ceiling(sqrt(p+1))
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nrow_plot, ncol_plot),
      mar = c(3.2, 3.5, 2.2, 0.8), oma = c(0, 0, 2.5, 0),
      mgp = c(2, 0.7, 0))
  
  # Plot each parameter
  x <- seq_len(M)
  for (j in seq_len(p)) {
    y <- theta_mat[, j]
    # Clamp plotting range to [0,1] as requested
    plot(x, y, type = "l", lwd = lwd_trace, col = col_trace,
         #ylim = c(0, 1),
         xlab = "Iteration", ylab = "",
         main =bquote(theta[.(j)]))
    abline(h = theta_true[j], col = col_true, lty = lty_true, lwd = 2)
    # Optional: show minor grid
    grid(col = adjustcolor("grey80", 0.8))
  }
  if (add_I){
  plot(simI,type="l",col=col_trace,lwd=lwd_trace,xlab="Time",ylab="Incidence")
  lines(trueI,type="l",lty=lty_true,col=col_true,lwd=2)
  grid(col = adjustcolor("grey80", 0.8))
  }
  mtext(main, outer = TRUE, cex = 1.1, font = 2)
}


plot_latent_traces <- function(x_mat, x_true,dimensions=2:8,
                              main = "Parameter trace plots",
                              col_trace = "#1f77b4", col_true = "#d62728",
                              lwd_trace = 1.5, lty_true = 2,add_I=T) {
  
  p=dim(x_mat)[3]-1
  m=length(dimensions)
  # Layout: try to make it roughly square
  ncol_plot <- 2#ceiling(sqrt(p+1))
  nrow_plot <- ceiling((m+1)/2)#ceiling(sqrt(p+1))
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nrow_plot, ncol_plot),
      mar = c(3.2, 3.5, 2.2, 0.8), oma = c(0, 0, 2.5, 0),
      mgp = c(2, 0.7, 0))
  
  # Plot each parameter
  
  for (j in seq_len(m)) {
    y <- x_mat[,dimensions[j],p+1]
    # Clamp plotting range to [0,1] as requested
    plot(y, type = "l", lwd = lwd_trace, col = col_trace,
         xlab = "Iteration", ylab = "",
         main = bquote(X[.(dimensions[j])]))
    abline(h = x_true[dimensions[j]], col = col_true, lty = lty_true, lwd = 2)
    # Optional: show minor grid
    grid(col = adjustcolor("grey80", 0.8))
  }
  mtext(main, outer = TRUE, cex = 1.1, font = 2)
}

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

  esjd <- function(arr) {
 
    # arr: array of dimension [it, nT, d]
    
    dims <- dim(arr)
    it <- dims[1]
    nT <- dims[2]
    d  <- dims[3]
    
    # storage for ESJD per nT
    esjd <- numeric(nT)
    
    for (j in 1:nT) {
      # extract trajectory for this nT: matrix [it, d]
      traj <- arr[, j,1:d ]
      
      # compute differences between consecutive iterations
      diffs <- traj[2:it,1:d , drop = FALSE] - traj[1:(it-1), 1:d, drop = FALSE]
      
      # squared Euclidean distance at each step
      sq_jump <- rowSums(diffs^2)
      
      # average over iterations
      esjd[j] <- mean(sq_jump)
    }
    
    return(esjd)
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