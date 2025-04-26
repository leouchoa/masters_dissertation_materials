#' Log-likelihood of the Bivariate (Wittle-)Matern Model.
#'
#' @param obs_matrix A n x 2 numeric matrix of data. We assume that the user has already log-transformed the data.
#' @param coords_matrix A n x 2 numeric matrix of coordinates for the data.
#' @param theta0 Initial values for the optimizer. A 4 x 1 numeric vector, the order of parameters is: theta = c(sigma_1^2, sigma_2^2, a, rho).
#' @param nus A 2 x 1 numeric vector with the smoothness coefficients for the Mat\'{e}rn covariance function, c(nu1, nu2).
#'
#' @return \code{fit_biwm} returns an object of \code{\link{class}} "biwm" containing at least the following components
#'   \item{mu}{Expected value estimate}
#'   \item{theta}{TODO}
#'   \item{se.mu}{TODO}
#'   \item{se.theta}{TODO}
#'   \item{likelihood}{TODO}
#'
#' @export
#'
#' @examples
#' set.seed(3)
#' n <- 40
#' coords <- matrix(runif(2*n), ncol = 2)
#' temp <- rnorm(2*n)
#' S <- sigma_assembler_biwm(sigmas = c(1, 1), a = 2, rho = 0.5, nus = c(0.5, 0.5), coords_matrix = coords)
#' log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)
#'
#' informed_test <- fit_biwm(log_cd, coords, c(1, 1, 2, .5), c(0.5, 0.5))
#' generic_test <- fit_biwm(log_cd, coords, c(.5, .5, 4, .6), c(0.5, 0.5))
#'
fit_biwm <- function(obs_matrix,
                     coords_matrix,
                     theta0,
                     nus, verbosity = 0, ...){

  if(length(nus) != 2){
    stop("Smoothness parameter vector must be of length 2", call. = FALSE)
  }

  if(length(theta0) != 4){
    stop("All shape paramaters must be the same", call. = FALSE)
  }

  if(abs(theta0[4]) > sqrt(nus[1]*nus[2])/mean(nus)){
    stop("Rho inserted does not define a valid covariance matrix", call. = FALSE)
  }

  # Average smoothness for third entry
  nus[3] <- mean(nus[1:2])

  mu <- colMeans(obs_matrix) # maximum profile-likelihood, initial value

  # LLike_biwm(theta0, nus = nus,  mu = mu, coords_matrix = coords_matrix, obs_matrix = obs_matrix)
  # LLike_biwm_reduced_grad(theta0, nus = nus,  mu = mu, coords_matrix = coords_matrix, obs_matrix = obs_matrix)
  mle <- optim(theta0,
               fn = LLike_biwm,
               gr = block_LLike_biwm_grad,
               method = "L-BFGS-B",
               lower = c(0.001, 0.001, 0.001, -min(sqrt(nus[1]*nus[2])/mean(nus), 0.999)),
               upper = c(Inf, Inf, Inf, min(sqrt(nus[1]*nus[2])/mean(nus), 0.999)),
               control = list(fnscale = -1, trace = verbosity),
               nus = nus,  mu = mu, # We always use Z = vec(\mathbf{Y} - boldsymbol\mu), maybe center it outside function?
               coords_matrix = coords_matrix,
               obs_matrix = obs_matrix)

  # Update mu with generalized least squares, maybe update mle

  ret <- list(mu = mu,
              theta = mle$par,
              likelihood = mle$value)

  class(ret) <- "biwm"

  return(ret)

}
