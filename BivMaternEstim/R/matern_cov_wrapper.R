# Matern covariance wraper for coordinate matrices
#
# Future work: generalize to cov_wrapper instead of matern_cov, see fields options.
matern_cov_wrapper <- function(coords_dist, a, nu){

  corr <- as.matrix(fields::Matern(d = coords_dist,
                                   alpha = a,
                                   nu = nu))
  return(corr + diag(min(nrow(corr), ncol(corr))))

}
