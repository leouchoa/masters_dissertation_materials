#' Log-likelihood of the Bivariate (Wittle-)Matern Model.
#'
#' @param obs_vec A n x 2 numeric matrix of data. We assume that the user has already log-transformed the data.
#' @param coords_matrix A n x 2 numeric matrix of coordinates for the data.
#' @param theta0 Initial values for the optimizer. A 4 x 1 numeric vector, the order of parameters is: theta = c(sigma_1^2, sigma_2^2, a, rho).
#' @param nus A 2 x 1 numeric vector with the smoothness coefficients for the Mat\'{e}rn covariance function, c(nu1, nu2).
#'
#' @return \code{fit_biwm} returns an object of \code{\link{class}} "biwm" containing at least the following components
#'   \item{mu}{Expected value estimate}
#'   \item{theta}{TODO}
#'   \item{se.mu}{TODO}
#'   \item{se.theta}{TODO}
#'   \item{AIC}{TODO}
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 40
#' coords <- matrix(runif(2*n), ncol = 2)
fit_biwm <- function(obs_vec,
                     coords_matrix,
                     theta0,
                     nus){

  mu <- colMeans(obs_vec) # maximum profile-likelihood, initial value

  ret <- list(mu = mu,
              theta = theta)

  class(ret) <- "biwm"

  return(ret)

}
