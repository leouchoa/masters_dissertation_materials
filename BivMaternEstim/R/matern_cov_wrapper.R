# Matern covariance wraper for coordinate matrices
#
# Future work: generalize to cov_wrapper instead of matern_cov, see fields options.
matern_cov_wrapper <- function(coords_dist, a, nu){

  as.matrix(
    fields::Matern(d = coords_matrix,
                   alpha = a,
                   nu = nu)
  )
}
