#' Compute the Conditional Multivariate Normal Covariance Matrix
#'
#'
#' @examples
#'
#'set.seed(1234)
#'loc_obs <- matrix(runif(2*10),ncol =2)
#'loc_new <- matrix(runif(2*3),ncol =2)
#'
#'sigmas <- c(1,2)
#'a <- 0.3
#'rho <- 0.8
#'nus <- c(1,1)
#'nug_vec <- c(0,0)
#'
#'
#'std_err_mat <- conditional_cov_mat(loc_new,loc_obs,sigmas, a, rho, nus, nug_vec)
#'
#'
conditional_cov_mat <- function(loc_new,loc_obs,sigmas, a, rho, nus, nug_vec){


  d_loc_obs <- dist(loc_obs)

  Sigma_hat_11 <-
    sigma_assembler_biwm(sigmas = sigmas,
                         a = a,
                         rho = rho,
                         nus = nus,
                         coords_matrix = loc_new,
                         nug_vec = nug_vec,
                         combined = TRUE)

  Sigma_hat_22 <-
    sigma_assembler_biwm(sigmas = sigmas,
                         a = a,
                         rho = rho,
                         nus = nus,
                         coords_matrix = loc_obs,
                         nug_vec = nug_vec,
                         combined = TRUE)

  Sigma_hat_12 <- cov_with_new_loc(loc_new = loc_new,
                   loc_obs = loc_obs,
                   sigmas = sigmas,
                   a = a,
                   rho = rho,
                   nus = nus,
                   nug_vec = nug_vec)



  return(
    Sigma_hat_11 - Sigma_hat_12 %*% solve(Sigma_hat_22,t(Sigma_hat_12))
  )

}
# source("R/matern_cov_wrapper.R")
# source("R/cov_with_new_loc.R")
# source("R/matern_cov_wrapper_krig.R")
