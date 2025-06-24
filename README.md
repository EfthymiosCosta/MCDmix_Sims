# MCDmix
## Extending the MCD estimator to mixed continuous-ordinal data

This repository includes functions that implement the extension of the Minimum Covariance Determinant (MCD) estimator to mixed-type data. The repository is organised as follows:

* `data`: Directory that includes the London Airbnb listings data file in `.csv` format. The data set can also be found in [Kaggle](https://www.kaggle.com/datasets/thedevastator/airbnb-prices-in-european-cities/data)
* `res`: Results directory. Includes the results files obtaied by running the scripts on the main repositoy page in `.RDS` format:
* `src`: Main functions for implementing MCDmix:
    * `find_min_lambda.R`: Bisection algorithm for finding the regularisation strength `lambda` that ensures the covariance matrix is positive definite with trace of at most `max_kappa` up to precision `tol`.
    * `lspace_proj_ord.R`: Function projecting ordinal variables to latent Gaussian space using the expectation of a truncated Gaussian with support defined by the manifest ordinal levels.
    * `mcd_fun_ord.R`: Main implementation of MCDmix for mixed continuous-ordinal data. Functions for computing the critical value for outlier identification, the mean squared error and the Kullback-Leibler divergence between two covariance matrices are also included.
* `MCDmix_Airbnb.R`: Script for running MCDmix on the London Airbnb data set.
* `MCDord_clustcont.R`: Script for running simulation study with outliers defined by cluster contamination.
* `MCDord_infreq_1cat.R`: Script for running simulation study with one highly infrequent ordinal level.
* `MCDord_mixedcorr.R`: Script for running simulation study with mixed correlation outliers defined between continuous and an ordinal variable.
    * Note that in the `res` directory, there are 4 output files corresponding to 2 up to 5 levels for the ordinal level that is defined as a linear combination of the continuous features. The code is identical for all cases, except that the number of levels changes.
* `MCDord_shift.R`: Script for running simulation study with outliers defined by shift contamination.
