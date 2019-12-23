# Matern covariance wraper for coordinate matrices
matern_cov_wrapper <- function(coords_matrix, a, nu){

  matrix(
    fields::Matern(d = as.vector(as.matrix(dist(coords_matrix))),
                   alpha = a,
                   nu = nu),
    ncol = nrow(coords_matrix),
    nrow = nrow(coords_matrix)
  )
}
