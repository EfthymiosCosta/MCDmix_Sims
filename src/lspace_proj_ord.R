lspace_proj_ord <- function(data, cov_mat, cont_cols, cuts, lvls) {
  
  n <- nrow(data)
  n_cont <- length(cont_cols)
  p <- ncol(cov_mat)
  
  # Precompute constants
  inv_cov <- solve(cov_mat[1:n_cont, 1:n_cont])
  cm <- colMeans(data[, cont_cols])
  
  # Precompute the matrix used to compute new_vec
  M <- cov_mat[(n_cont + 1):(p), 1:n_cont] %*% inv_cov
  
  # Precompute the constant new covariance and standard deviations
  new_cov <- cov_mat[(n_cont + 1):(p), (n_cont + 1):(p)] -
    cov_mat[(n_cont + 1):(p), 1:n_cont] %*%
    inv_cov %*% cov_mat[1:n_cont, (n_cont + 1):(p)]
  sd_vals <- sqrt(diag(new_cov))
  
  # Initialise result
  df_res <- matrix(NA, nrow = n, ncol = p-n_cont)
  
  # Loop over rows
  for (i in 1:n) {
    # Compute continuous deviation
    diff <- as.numeric(data[i, cont_cols] - cm)
    new_vec <- as.numeric(M %*% diff)
    
    # Thresholds
    x <- as.numeric(data[i, (n_cont + 1):(p)])
    supports_lower <- c()
    supports_upper <- c()
    for (j in 1:length(x)){
      supports_lower <- c(supports_lower,
                          ifelse(x[j] == 1, -Inf, cuts[[j]][x[j] - 1]))
      supports_upper <- c(supports_upper,
                          ifelse(x[j] == lvls[j], Inf, cuts[[j]][x[j]]))
    }
    
    # Project
    vals <- truncnorm::etruncnorm(a = supports_lower,
                                  b = supports_upper,
                                  mean = new_vec,
                                  sd = sd_vals)
    # Assign
    df_res[i, ] <- vals
  }
  
  return(as.data.frame(df_res))
}
