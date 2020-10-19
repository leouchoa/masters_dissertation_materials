#' Obtain Prediction Standard Error Estimate
#'
#' @param biwm_fit An object of biwm class. Must be the result of a `fit_biwm` command.
#'
#' @param loc_new A n x 2 numeric matrix of coordinates for the unobserved data.
#'
#' @param loc_obs A n x 2 numeric matrix of coordinates for the estimation data.
#'
#' @param nus A 2 x 1 numeric vector with the smoothness coefficients for the Mat\'{e}rn covariance function, c(nu1, nu2).
#'
#' #' @param nug_vec A 2 x 1 numeric vector with nugget effects Mat\'{e}rn covariance function, c(nu1, nu2).
#'
#'
#' @export
#' @examples
#' set.seed(3)
#' n <- 40
#' coords <- matrix(runif(2*n), ncol = 2)
#' temp <- rnorm(2*n)
#' nug_vec <- c(0,0)
#' nus <- c(0.5,0.5)
#' S <- sigma_assembler_biwm(sigmas = c(1, 1), a = 2, rho = 0.5,
#' nus = nus, coords_matrix = coords, nug_vec = nug_vec, combined = TRUE)
#' log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)
#'
#' generic_test <- fit_biwm(log_cd, coords, c(.5, .5, 4, .6), nus = nus,nug_vec = nug_vec)
#'
#' loc_new <- matrix(runif(2*20),ncol =2)
#'
#' pred_std_err_mat <- get_pred_std_err(loc_new = loc_new, loc_obs = coords,biwm_fit = generic_test, nus = nus, nug_vec = nug_vec)
#'
#'
get_pred_std_err <- function(loc_new,loc_obs,biwm_fit,nus,nug_vec){

  cond_cov_mat <- conditional_cov_mat(loc_new = loc_new,
                                      loc_obs = loc_obs,
                                      sigmas = biwm_fit$theta[1:2],
                                      a = biwm_fit$theta[3],
                                      rho = biwm_fit$theta[4],
                                      nus = nus,
                                      nug_vec = nug_vec)

  aux_df <-
    alr_inv(
      cbind(
        diag(cond_cov_mat)[1:(nrow(loc_new))],
        diag(cond_cov_mat)[(nrow(loc_new) + 1):(2*nrow(loc_new))]
      )
    )

  return(
    cbind(setNames(as.data.frame(loc_new),c("coord_x","coord_y")),
          comp_01_pred_std_err = aux_df[,1],
          comp_02_pred_std_err = aux_df[,2],
          comp_03_pred_std_err = aux_df[,3]
          )
  )
}

