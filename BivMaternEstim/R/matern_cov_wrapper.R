#Fields implementation at cara dura

Matern <- function (d, alpha = 1/range, nu = 0.5){
  
  if (any(d < 0)) 
    stop("distance argument must be nonnegative")
  
  d <- d * alpha
  d[d == 0] <- 1e-10
  con <- (2^(nu - 1)) * gamma(nu)
  con <- 1/con
  return(con * (d^nu) * besselK(d, nu))
}

# Matern covariance wraper for coordinate matrices

matern_cov_wrapper <- function(coords_matrix, a, nu){
  
  matrix(
    Matern(d = as.vector(as.matrix(dist(coords_matrix))),
           alpha = a, 
           nu = nu),
    ncol = nrow(coords_matrix),
    nrow = nrow(coords_matrix)
  )
}
