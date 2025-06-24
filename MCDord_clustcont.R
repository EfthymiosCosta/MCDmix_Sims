source('src/mcd_fun_ord.R')
source('src/find_min_lambda.R')
source('src/lspace_proj_ord.R')
library(clusterGeneration)
# categorized numerical variable function
intv <- function(vec, class) {
  nbase <- (1:(class-1))/class
  nq <- numeric(length(nbase))
  for (i in 1:length(nq)) {
    nq[i] <- quantile(vec, nbase[i])
  }
  res <- c(min(vec), nq, max(vec)) 
  res[1] <- res[1]-1
  for (i in 2:length(res)){
    if (res[i-1]==res[i]){
      res[i] <- res[i]+2e-15
    }
  }
  return(res)
}

n <- 500
pC <- 4
seeds <- c(1:100)
eps_vec <- c(0.05, 0.10)
h_vec <- c(0.75, 0.8, 0.9)
lvls_vec <- c(3, 4, 2)
n_sim <- 7*length(seeds)
clust_cont_outs_ord <- vector("list", n_sim)
clust_cont_outs_ord <- lapply(clust_cont_outs_ord, function(x) list())
count <- 1

for (seed_num in seeds){
  for (eps_val in eps_vec){
    set.seed(seed_num)
    num_outs <- round(n*eps_val)
    p <- pC + length(lvls_vec)
    # Generate random covariance matrix for continuous variables
    rand_sigma <- clusterGeneration::genPositiveDefMat(dim = p,
                                                       covMethod = "unifcorrmat",
                                                       alphad = 1,
                                                       rangeVar = c(1, 1))$Sigma
    # Generate normal (non-outlying) data
    random_dt <- mvtnorm::rmvnorm(n = 500 - num_outs,
                                  mean = rep(0, p),
                                  sigma = rand_sigma)
    random_dt <- data.frame(random_dt)
    # Now take eigenvector corresponding to smallest eigenvalue and rescale
    min_evec <- eigen(rand_sigma)$vectors[, p]
    const <- c(t(min_evec) %*% solve(rand_sigma) %*% min_evec)
    min_evec <- sqrt(p)*min_evec/sqrt(const)
    
    # Generate these outliers setting k = 50
    random_dt_2 <- mvtnorm::rmvnorm(n = num_outs,
                                    mean = c(50*min_evec),
                                    sigma = 0.05^2*diag(p))
    random_dt_2 <- data.frame(random_dt_2)
    random_dt <- rbind(random_dt, random_dt_2)
    # Set outlier indices
    out_indices <- c((n - num_outs + 1):n)
    # Convert to ordinal
    for (j in c((pC+1):p)){
      random_dt[, j] <- as.factor(cut(random_dt[, j],
                                      intv(random_dt[ , j], (lvls_vec[j-pC])),
                                      labels = (1:(lvls_vec[j-pC]))))
    }
    for (h_val in h_vec){
      if (eps_val + h_val > 1) next
      # Run MCD
      MCDmix <- mcd_mixed_ord(dt = random_dt,
                              cont_cols = c(1:pC),
                              reps = 25,
                              h = h_val,
                              max_kappa = 50)
      # Compute critical value
      crit_val <- outlier_detect_critical_val(alpha = 0.01,
                                              n = nrow(random_dt),
                                              dfs = ncol(random_dt))
      # Add parameters in output list
      clust_cont_outs_ord[[count]]$h <- h_val
      clust_cont_outs_ord[[count]]$epsilon <- eps_val
      # Calculate proportion of outliers in h-subset
      clust_cont_outs_ord[[count]]$prop_in_h <- length(intersect(MCDmix$Subset, out_indices))/num_outs
      # Calculate proportion of outliers detected by algorithm
      outs_detected <- which(MCDmix$Mahalanobis > crit_val)
      clust_cont_outs_ord[[count]]$prop_true_outs_detected <- length(intersect(outs_detected, out_indices))/num_outs
      # Calculate proportion of regular observations marked as outliers
      clust_cont_outs_ord[[count]]$prop_false_outs_detected <- length(setdiff(outs_detected, out_indices))/length(outs_detected)
      # Number of outliers among the n*(1-epsilon) observations with the lowest distance
      clust_cont_outs_ord[[count]]$num_smallest_dist <- length(intersect(out_indices, order(MCDmix$Mahalanobis)[1:(n-num_outs)]))
      # Compute KL and MSE
      clust_cont_outs_ord[[count]]$MSE <- mse_mat(rand_sigma, MCDmix$Covariance)
      clust_cont_outs_ord[[count]]$KL <- kl_mat(rand_sigma, MCDmix$Covariance)
      
      saveRDS(clust_cont_outs_ord, file = 'clust_cont_outs_ord.RDS')
      cat('Count', count,'/', n_sim, 'complete. :) \n')
      count <- count + 1
    }
  }
}
