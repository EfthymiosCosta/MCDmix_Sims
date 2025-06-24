find_min_lambda <- function(cov_mat,
                            reg_mat,
                            lambda_low,
                            lambda_high,
                            max_kappa,
                            tol) {
  # Ensure symmetry for initial matrices
  cov_mat <- as.matrix(Matrix::forceSymmetric(cov_mat))
  reg_mat <- as.matrix(Matrix::forceSymmetric(reg_mat))
  
  # Bisection loop
  while ((lambda_high - lambda_low) > tol){
    lambda_mid <- (lambda_low + lambda_high) / 2
    cov_mat_reg <- (1 - lambda_mid) * cov_mat + lambda_mid * reg_mat
    cov_mat_reg <- as.matrix(Matrix::forceSymmetric(cov_mat_reg))
    
    if (matrixcalc::is.positive.definite(cov_mat_reg) && kappa(cov_mat_reg, exact = TRUE) <= max_kappa){
      # If positive definite and kappa <= max_kappa, try a smaller lambda
      lambda_high <- lambda_mid
    } else {
      # If not positive definite or kappa > max_kappa, increase lambda
      lambda_low <- lambda_mid
    }
  }
  
  # Return smallest lambda
  return(list(lambda_high, kappa(cov_mat_reg, exact = TRUE)))
}
