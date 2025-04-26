# Matern covariance wraper for coordinate matrices
#
# Future work: generalize to cov_wrapper instead of matern_cov, see fields options.
#' @importFrom fields Matern
matern_cov_wrapper_krig <- function(coords_dist, a, nu){
  
  corr <- as.matrix(Matern(d = coords_dist,
                           alpha = a,
                           nu = nu))
}
