#' Log-likelihood of the Bivariate (Wittle-)Matern Model.
#'
#' @param theta A 4 x 1 numeric vector, the order of parameters is: theta = c(sigma_1^2, sigma_2^2, a, rho).
#' @param nus A 2 x 1 numeric vector with the smoothness coefficients for the Mat\'{e}rn covariance function, c(nu1, nu2).
#' @param mu A 2 x 1 numeric vector of means c(mu1, mu2).
#' @param coords_matrix A n x 2 numeric matrix of coordinates for the data.
#' @param obs_vec A n x 2 numeric matrix of data. We assume that the user has already log-transformed the data.
#'
#' @return log-likelihood
#'
#' @export
LLike_biwm <- function(theta,
                       nus,
                       mu,
                       coords_matrix,
                       obs_vec){

  sigmas <- theta[1:2]
  as <- theta[3]
  rho <- theta[4]

  if(length(nus) != 2){
    stop("Smoothness parameter vector must be of length 2", call. = FALSE)
  }

  if(length(unique(as)) != 1){
    stop("All shape paramaters must be the same", call. = FALSE)
  }

  if(abs(rho) > sqrt(nus[1]*nus[2])/mean(nus)){
    stop("Rho inserted does not define a valid covariance matrix", call. = FALSE)
  }

  # Average smoothness for third entry
  nus[3] <- mean(nus[1:2])

  if(length(as) == 1){
    as <- rep(as, 3)
  }

  autocov_matrix <- sigma_assembler_biwm(sigmas = sigmas,
                                         as = as,
                                         rho = rho,
                                         nus = nus,
                                         coords_matrix = coords_matrix)

  as.double(
    -1/2 *
      (
        determinant(autocov_matrix, logarithm = TRUE) +
          crossprod(obs_vec - mu, solve(autocov_matrix, obs_vec - mu)) +
          length(obs_vec) * log(2*pi)
      )
  )

}
