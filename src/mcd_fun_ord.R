source('find_min_lambda.R')
source('lspace_proj_ord.R')
source('lspace_proj_ord_nocont.R')
library(polycor)
library(caret)
library(propagate)
library(matrixStats)
library(matrixcalc)
library(truncnorm)
library(robustbase)

outlier_detect_critical_val <- function(alpha, n, dfs){
  alpha_n <- (1 - alpha)^(1/n)
  critical_val <- qchisq(p = alpha_n, df = dfs)
  return(critical_val)
}

mse_mat <- function(mat1, mat2){
  p <- nrow(mat1)
  return(sum((mat1 - mat2)^2)/p^2)
}

kl_mat <- function(mat1, mat2){
  # Mat2 must be the estimated matrix
  p <- nrow(mat1)
  Shat_Sinv <- mat2 %*% solve(mat1)
  return(sum(diag(Shat_Sinv)) - log(det(Shat_Sinv)) - p)
}

mcd_mixed_ord <- function(dt, cont_cols, reps = 10,
                          h = 0.75, max_kappa = 50,
                          max_iter = 50){
  dt <- as.data.frame(dt)
  # Scaling using median and MAD
  for (i in cont_cols){
    dt[, i] <- (dt[, i] - median(dt[, i])) / mad(dt[, i])
  }
  cat_cols <- setdiff(c(1:ncol(dt)), cont_cols)
  for (j in cat_cols){
    dt[, j] <- as.numeric(dt[, j])
  }
  cat_vars <- rep(1, length(cat_cols))
  num_lvls <- as.numeric(apply(data.frame(dt[, cat_cols]), 2,
                         FUN = function(x) length(unique(x))))
  # Vector to store row cuts for contingency table for tetrachoric correlation
  cuts <- vector(mode = "list", length = length(cat_cols))
  full_corr_mat <- matrix(NA, nrow = ncol(dt), ncol = ncol(dt))
  full_corr_mat[c(1:length(cont_cols)), c(1:length(cont_cols))] <- cor(as.matrix(dt[, cont_cols]))
  diag(full_corr_mat) <- 1
  list_count <- 1
  # Polychoric correlations between binary features
  if (max(cont_cols) == (ncol(full_corr_mat) - 1)){
    poly_obj <- polychor(dt[, ncol(full_corr_mat)],
                         dt[, ncol(full_corr_mat)],
                         ML = FALSE,
                         thresholds = TRUE,
                         maxcor = 1)
    cuts[[list_count]] <- c(cuts[[list_count]], poly_obj$row.cuts)
  } else {
    for (i in (length(cont_cols)+1):(ncol(full_corr_mat))){
      for (j in setdiff((length(cont_cols)+1):(ncol(full_corr_mat)), i)){
        # Need to check for constant columns
        if (length(unique(dt[, i])) == 1 | length(unique(dt[, j])) == 1){
          full_corr_mat[i, j] <- 0
        } else {
          poly_obj <- polychor(dt[, i],
                               dt[, j],
                               ML = FALSE,
                               thresholds = TRUE,
                               maxcor = 1)
          full_corr_mat[i, j] <- poly_obj$rho
        }
      }
      cuts[[list_count]] <- c(cuts[[list_count]], poly_obj$row.cuts)
      list_count <- list_count + 1
    }
  }
  
  # Polyserial correlations between continuous and binary features
  for (i in (length(cont_cols)+1):(ncol(full_corr_mat))){
    for (j in 1:length(cont_cols)){
      # Need to check for constant columns
      if (length(unique(dt[, i])) == 1 | length(unique(dt[, j])) == 1){
        full_corr_mat[i, j] <- 0
        full_corr_mat[j, i] <- 0
      } else {
        full_corr_mat[i, j] <- polyserial(y = dt[, i],
                                          x = dt[, j],
                                          ML = FALSE,
                                          thresholds = TRUE,
                                          maxcor = 1)$rho
        full_corr_mat[j, i] <- full_corr_mat[i, j]
      }
    }
  }
  
  # Convert to covariance matrix
  full_cov_mat <- propagate::cor2cov(full_corr_mat,
                                     var = c(matrixStats::colVars(as.matrix(dt[,cont_cols])),
                                             cat_vars))
  # Convert binary features using the tetrachoric correlation cut-off points
  # Use the mean of the truncated normal distributions
  dt_latent_space <- dt
  
  # Diagonal matrix with MADs and Variances for regularisation
  ### ALTERNATIVE REGULARISATION MATRIX FOR DUMMIES:
  ### Take 1 + p(1-p)
  reg_mat <- diag(c(colMads(as.matrix(dt[,c(1:length(cont_cols))]))^2,
                    cat_vars))
  # Consistency factor for covariance matrix
  consistency_factor <- h/(pchisq(qchisq(h, df = ncol(dt)), df = ncol(dt)+2))
  
  best_sub <- list('Subset'= c(NA), 'Determinant' = Inf, 'Covariance' = matrix(NA, nrow = ncol(dt), ncol = ncol(dt)))
  for (rep in 1:reps){
    cat('Rep', rep,'\n')
    
    # Start with initial subset
    set.seed(rep)
    init_sub <- sample(c(1:nrow(dt)), length(cont_cols) + 1)
    
    # Compute initial Mahalanobis distances, then apply C-step
    Sigma_init <- cov(as.matrix(dt[init_sub, cont_cols]))
    mu_init <- colMeans(as.matrix(dt[init_sub, cont_cols]))
    md_init <- mahalanobis(as.matrix(dt[, cont_cols]), center = mu_init, cov = Sigma_init)
    #md_init[which(dt$V3=='3')] <- Inf
    subset <- order(md_init)[1:ceiling(h*nrow(dt))]
    count <- 1
    iter <- 1
    old_det <- Inf
    # Switch variable
    aux <- TRUE
    while (aux){
      # Initialise correlation matrix
      corr_mat <- matrix(NA, nrow = ncol(dt), ncol = ncol(dt))
      # Compute correlation matrix elements for continuous variables (on subset)
      corr_mat[c(1:length(cont_cols)), c(1:length(cont_cols))] <- cor(as.matrix(dt[subset, cont_cols]))
      # Set diagonal elements to 1
      diag(corr_mat) <- 1
      
      # Polychoric correlations between binary features
      for (i in (length(cont_cols)+1):(ncol(corr_mat))){
        for (j in setdiff((length(cont_cols)+1):(ncol(corr_mat)), i)){
          # Need to check for constant columns
          if (length(unique(dt[subset, i])) == 1 | length(unique(dt[subset, j])) == 1){
            corr_mat[i, j] <- 0
          } else {
            tbl_sub <- table(dt[subset, i], dt[subset, j])
            # Apply correction
            #if (sum(tbl_sub == 0) == 1){
            #  tbl_sub[which(tbl_sub==0)] <- 0.5
            #}
            #if (sum(tbl_sub == 0) == 2){
            #  if (tbl_sub[1] == 0){
            #    corr_mat[i, j] <- -1
            #  } else if (tbl_sub[2] == 0){
            #    corr_mat[i, j] <- 1
            #  }
            #} else {
              corr_mat[i, j] <- polychor(tbl_sub,
                                         ML = FALSE,
                                         thresholds = TRUE,
                                         maxcor = 1)$rho
            #}
          }
        }
      }
      
      # Polyserial correlations between continuous and binary features
      for (i in (length(cont_cols)+1):(ncol(corr_mat))){
        for (j in 1:length(cont_cols)){
          # Need to check for constant columns
          if (length(unique(dt[subset, i])) == 1 | length(unique(dt[subset, j])) == 1){
            corr_mat[i, j] <- 0
            corr_mat[j, i] <- 0
          } else {
            corr_mat[i, j] <- polyserial(y = dt[subset, i],
                                         x = dt[subset, j],
                                         ML = FALSE,
                                         thresholds = TRUE,
                                         maxcor = 1)$rho
            corr_mat[j, i] <- corr_mat[i, j]
          }
        }
      }
      
      # Convert to covariance matrix
      cov_mat <- propagate::cor2cov(corr_mat, var = c(matrixStats::colVars(as.matrix(dt[subset, cont_cols])),
                                                      cat_vars))
      # Force symmetry to avoid numerical issues
      cov_mat <- as.matrix(Matrix::forceSymmetric(cov_mat))
      # Check if matrix is positive semi-definite
      # If not, need to regularise
      lambda_search_aux <- FALSE
      kappa_val <- kappa(cov_mat, exact = TRUE)
      if (!matrixcalc::is.positive.definite(cov_mat*consistency_factor) || kappa_val > max_kappa){
        lambda_search <- find_min_lambda(cov_mat*consistency_factor,
                                         reg_mat,
                                         lambda_low = 0,
                                         lambda_high = 1,
                                         max_kappa = max_kappa,
                                         tol = 1e-4)
        lambda <- lambda_search[[1]]
        kappa <- lambda_search[[2]]
        cov_mat_reg <- (1-lambda)*consistency_factor*cov_mat + lambda*reg_mat
        lambda_search_aux <- TRUE
        #cat('Lambda =', lambda, '\n')
        #cat('Kappa =', kappa, '\n')
      } else {
        cov_mat_reg <- cov_mat*consistency_factor
      }
      cov_mat_reg <- as.matrix(Matrix::forceSymmetric(cov_mat_reg))
      # Project binary features to latent space with updated covariance matrix
      dt_latent_space[, c((length(cont_cols)+1):(ncol(dt)))] <- lspace_proj_ord(dt, cov_mat_reg, cont_cols, cuts, num_lvls,
                                                                                cm = colMeans(dt[subset, c(1:length(cont_cols))]))
      # Compute Mahalanobis distances
      mu <- c(colMeans(as.matrix(dt_latent_space[subset, cont_cols])), rep(0, ncol(dt)-length(cont_cols)))
      mds <- mahalanobis(dt_latent_space, center = mu, cov = cov_mat_reg)
      # Use the C-step
      new_subset <- order(mds)[1:ceiling(h*nrow(dt))]
      # Check if subset remains unchanged
      if (all(new_subset %in% subset | iter == max_iter)){
        old_det <- det(cov_mat_reg)
        old_cov_mat_reg <- cov_mat_reg
        aux <- FALSE
      } else {
        old_subset <- subset
        subset <- new_subset
        count <- count + 1
        old_det <- det(cov_mat_reg)
        old_cov_mat_reg <- cov_mat_reg
      }
      iter <- iter + 1
    }
    if (old_det < best_sub[[2]]){
      best_sub[[1]] <- subset
      best_sub[[2]] <- old_det
      best_sub[[3]] <- old_cov_mat_reg
      colnames(best_sub[[3]]) <- colnames(dt)
      rownames(best_sub[[3]]) <- colnames(dt)
      cat('Rep:', rep, '\n')
      ifelse(lambda_search_aux,
             best_sub$lambda <- lambda_search[[1]],
             best_sub$lambda <- 0)
      #cat('ITERATION:', rep, '\n')
    }
  }
  best_sub$Cuts <- cuts
  dt_latent_space[, c((length(cont_cols)+1):(ncol(dt)))] <- lspace_proj_ord(dt, best_sub[[3]], cont_cols, cuts, num_lvls,
                                                                            cm = colMeans(dt[best_sub[[1]], c(1:length(cont_cols))]))
  best_sub$dtLatent <- dt_latent_space
  best_sub$kappa <- kappa(best_sub[[3]], exact = TRUE)
  best_sub$Mahalanobis <- mahalanobis(dt_latent_space,
                                      center = c(colMeans(as.matrix(dt_latent_space[best_sub[[1]], cont_cols])),
                                                 rep(0, ncol(dt_latent_space)-length(cont_cols))),
                                      cov = best_sub[[3]])
  return(best_sub)
}
