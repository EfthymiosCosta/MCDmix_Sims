generate_covmat_ALYZ <- function(p, CN, tol = 1e-6, max_iter = 100){
  # generate eigenvals
  lambda <- c(1, sort(runif(p - 2, min = 1, max = CN)), CN)
  
  # QR decomposition
  Y <- matrix(rnorm(p * p), nrow = p, ncol = p)
  U <- qr.Q(qr(Y))
  
  # initial covariance matrix
  Sigma <- U %*% diag(lambda) %*% t(U)
  
  # iterate until kappa <= tol
  iter <- 0
  converged <- FALSE
  
  while (iter < max_iter && !converged) {
    iter <- iter + 1
    # convert to correlation matrix
    D_inv_sqrt <- diag(1 / sqrt(diag(Sigma)))
    R <- D_inv_sqrt %*% Sigma %*% D_inv_sqrt
    
    # eigendecomposition
    eig <- eigen(R, symmetric = TRUE)
    lambda_R <- eig$values
    lambda_min <- min(lambda_R)
    
    # CN of R
    current_CN <- kappa(R, exact = TRUE)
    
    # convergence check
    if (abs(current_CN - CN) < tol) {
      converged <- TRUE
    } else {
      # set largest eigenval = CN * lambda_min
      lambda_R[which.max(lambda_R)] <- CN * lambda_min
      Lambda_new <- diag(lambda_R)
      # sigma update
      Sigma <- eig$vectors %*% Lambda_new %*% t(eig$vectors)
    }
  }
  
  if (!converged) warning("Failed to converge within max iterations.")
  return(Sigma)
}