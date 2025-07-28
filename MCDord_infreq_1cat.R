source('src/mcd_fun_ord.R')
source('src/find_min_lambda.R')
source('src/lspace_proj_ord.R')
source('src/ALYZ.R')
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
occurences <- c(1:5)
h_vec <- c(0.75, 0.8, 0.9)
lvls_list <- list(c(2), c(3), c(4))
n_sim <- length(h_vec)*length(seeds)*length(occurences)*length(lvls_list)
infreqs_list_1cat_ord <- vector("list", n_sim)
infreqs_list_1cat_ord <- lapply(infreqs_list_1cat_ord, function(x) list())
count <- 1

for (seed_num in seeds){
  for (lvls_vec in lvls_list){
    for (num_occurences in occurences){
      set.seed(seed_num)
      # Generate random covariance matrix for continuous variables
      rand_sigma <- generate_covmat_ALYZ(p, CN = 100)
      random_dt <- mvtnorm::rmvnorm(n = 500,
                                    mean = rep(0, (pC+1)),
                                    sigma = rand_sigma)
      random_dt <- data.frame(random_dt)
      # Maximum num_occurences values in last column
      max_obs <- order(random_dt[, (pC+1)], decreasing = TRUE)[1:num_occurences]
      # Convert to ordinal
      if (lvls_vec > 2){
        random_dt[, (pC+1)] <- as.factor(cut(random_dt[, (pC+1)],
                                             intv(random_dt[ , (pC+1)], (lvls_vec-1)),
                                             labels = (1:(lvls_vec-1))))
      } else {
        random_dt[, (pC+1)] <- 1
      }
      
      random_dt[, (pC+1)] <- as.numeric(random_dt[, (pC+1)])
      random_dt[max_obs, (pC+1)] <- lvls_vec
      random_dt[, (pC+1)] <- as.factor(random_dt[, (pC+1)])
      
      for (h_val in h_vec){
        cat('Seed:', seed_num, '\t h:', h_val, '\t Occurences:', num_occurences,
            '\t Levels:', lvls_vec, '\n')
        # Define the infrequent level observations - let's call them outliers
        out_indices <- max_obs
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
        # Add parameters
        infreqs_list_1cat_ord[[count]]$h <- h_val
        infreqs_list_1cat_ord[[count]]$n_infreq <- num_occurences
        infreqs_list_1cat_ord[[count]]$p <- ncol(random_dt)
        infreqs_list_1cat_ord[[count]]$levels <- lvls_vec
        # Calculate proportion of outliers in h-subset
        infreqs_list_1cat_ord[[count]]$prop_in_h <- length(intersect(MCDmix$Subset, out_indices))/num_occurences
        # Calculate proportion of overall outliers detected by algorithm
        outs_detected <- which(MCDmix$Mahalanobis > crit_val)
        infreqs_list_1cat_ord[[count]]$prop_true_outs_detected <- length(intersect(outs_detected, out_indices))/num_occurences
        # Store Mahalanobis distance & critical value
        infreqs_list_1cat_ord[[count]]$mahalanobis <- MCDmix$Mahalanobis[out_indices]
        infreqs_list_1cat_ord[[count]]$critical_value <- crit_val
        # Number of outliers among the n*(1-epsilon) observations with the lowest distance
        infreqs_list_1cat_ord[[count]]$num_smallest_dist <- rank(MCDmix$Mahalanobis)[out_indices]
        saveRDS(infreqs_list_1cat_ord, file = 'infreqs_list_1cat_ord.RDS')
        cat('Count', count,'/', n_sim, 'complete. :) \n')
        count <- count + 1
      }
    }
  }
}
