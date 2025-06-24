source('src/mcd_fun_ord.R')
airbnb_london <- read.csv('data/london_weekdays.csv')
# Remove Observation
airbnb_london <- airbnb_london[, -1]
# Also remove attr_index and rest_index and keep their normalised versions
airbnb_london <- airbnb_london[, -c(14, 16)]
# Change guest_satisfaction overall to lie in groups [0, 50], [51, 60], [61, 70], ..., [91, 100]
breaks <- c(0, 50, 60, 70, 80, 90, 100)
labels <- 1:6
airbnb_london[,10] <- cut(airbnb_london[, 10],
                          breaks = breaks,
                          labels = labels,
                          right = TRUE,  
                          include.lowest = TRUE)
# Convert all categorical to ordinal
cat_cols <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
cont_cols <- setdiff(c(1:ncol(airbnb_london)), cat_cols)
for (i in cat_cols){
  airbnb_london[, i] <- as.factor(airbnb_london[, i])
}
# Take log of RealSum & overall guest satisfaction score & rest_index_norm & attr_index_norm
airbnb_london$realSum <- log(airbnb_london$realSum)
airbnb_london$rest_index_norm <- log(airbnb_london$rest_index_norm)
airbnb_london$attr_index_norm <- log(airbnb_london$attr_index_norm)
airbnb_london$metro_dist <- log(airbnb_london$metro_dist)

airbnb_london_1 <- airbnb_london
### IMPORTANT: Current implementation requires the continuous variables to appear first
# Therefore we re-arrange
airbnb_london_1 <- cbind(airbnb_london_1[, cont_cols], airbnb_london_1[, setdiff(c(1:ncol(airbnb_london_1)), cont_cols)])
cont_cols <- c(1:length(cont_cols))

# Run MCD treating categorical variables as ordinal
MCDOrd_Airbnb_London <- mcd_mixed_ord(dt = airbnb_london_1,
                                      cont_cols = cont_cols,
                                      reps = 100,
                                      h = 0.75,
                                      max_kappa = 50,
                                      max_iter = 50)

# Critical value
crit_val_airbnb <- outlier_detect_critical_val(alpha = 0.05,
                                               n = nrow(airbnb_london_1),
                                               dfs = ncol(airbnb_london_1))
outs_detected_airbnb <- which(MCDOrd_Airbnb_London$Mahalanobis > crit_val_airbnb)