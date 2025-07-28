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
lvls_vec <- c(3, 4)
n_sim <- 7*length(seeds)
mixedcorr_outs_ord <- vector("list", n_sim)
mixedcorr_outs_ord <- lapply(mixedcorr_outs_ord, function(x) list())
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
    # Get projections/direction  along which we cut
    evec <- eigen(rand_sigma[c(1:pC), c(1:pC)])$vectors[, 1]
    projs <- random_dt[, c(1:pC)] %*% evec
    # Variance of new variable introduced
    var_new <- as.numeric(t(evec) %*% rand_sigma[c(1:pC), c(1:pC)] %*% evec)
    corrs_vec <- c()
    for (j in 1:p){
      corrs_vec <- c(corrs_vec, sum(evec * rand_sigma[c(1:pC), j])/sqrt(var_new))
    }
    #random_dt <- data.frame(random_dt)
    rand_sigma <- cbind(rand_sigma, corrs_vec)
    rand_sigma <- rbind(rand_sigma, c(corrs_vec, 1))
    # Now generate outliers
    # Get projections/direction  along which we cut
    additional_binary_var <- ifelse(projs > 0, 2, 1)
    random_dt <- cbind(random_dt, additional_binary_var)
    colnames(random_dt) <- rep("", ncol(random_dt))
    # Randomly select num_outs observations and switch the last binary variable
    # Ensure we have sufficiently many observations
    n_valid <- length(which(abs(projs) > 1.5))
    # Also ensure we have enough of each kind
    obs_pos <- which(projs > 1.5)
    n_pos <- length(obs_pos)
    obs_neg <- which(projs < -1.5)
    n_neg <- length(obs_neg)
    if (n_valid > num_outs & n_pos >= round(num_outs/2) & n_neg >= round(num_outs/2)){
      out_indices_pos <- sample(obs_pos, round(num_outs/2))
      out_indices_neg <- sample(obs_neg, num_outs - round(num_outs/2))
      mixedcorr_outs_ord[[count]]$thresh <- 1.5
    } else {
      obs_pos <- which(projs > 1)
      obs_neg <- which(projs < -1)
      out_indices_pos <- sample(obs_pos, round(num_outs/2))
      out_indices_neg <- sample(obs_neg, num_outs - round(num_outs/2))
      mixedcorr_outs_ord[[count]]$thresh <- 1
    }
    # Define the outliers
    out_indices <- c(out_indices_pos, out_indices_neg)
    random_dt[out_indices_pos, ncol(random_dt)] <- 1
    random_dt[out_indices_neg, ncol(random_dt)] <- 2
    
    # Convert to ordinal
    for (j in c((pC+1):p)){
      random_dt[, j] <- as.factor(cut(random_dt[, j],
                                      intv(random_dt[ , j], (lvls_vec[j-pC])),
                                      labels = (1:(lvls_vec[j-pC]))))
    }
    random_dt <- data.frame(random_dt)
    
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
      mixedcorr_outs_ord[[count]]$h <- h_val
      mixedcorr_outs_ord[[count]]$epsilon <- eps_val
      # Calculate proportion of outliers in h-subset
      mixedcorr_outs_ord[[count]]$prop_in_h <- length(intersect(MCDmix$Subset, out_indices))/num_outs
      # Calculate proportion of outliers detected by algorithm
      outs_detected <- which(MCDmix$Mahalanobis > crit_val)
      mixedcorr_outs_ord[[count]]$prop_true_outs_detected <- length(intersect(outs_detected, out_indices))/num_outs
      # Calculate number of regular observations marked as outliers
      mixedcorr_outs_ord[[count]]$num_false_outs_detected <- length(setdiff(outs_detected, out_indices))
      # Number of outliers among the n*(1-epsilon) observations with the lowest distance
      mixedcorr_outs_ord[[count]]$num_smallest_dist <- length(intersect(out_indices, order(MCDmix$Mahalanobis)[1:(n-num_outs)]))
      # Compute MSE
      mixedcorr_outs_ord[[count]]$MSE <- mse_mat(rand_sigma, MCDmix$Covariance)
      
      saveRDS(mixedcorr_outs_ord, file = 'mixedcorr_outs_ord.RDS')
      cat('Count', count,'/', n_sim, 'complete. :) \n')
      count <- count + 1
    }
  }
}
