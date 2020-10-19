#' Compositional Krigging for the Bivariate (Wittle-)Matern Model.
#'
#' @param biwm_fit An object of biwm class. Must be the result of a `fit_biwm` command.
#'
#' @param krig_locations A n x 2 numeric matrix of coordinates for the unobserved data.
#'
#' @param fit_locations A n x 2 numeric matrix of coordinates for the estimation data.
#'
#' @param obs_matrix A n x 2 numeric matrix of data. We assume that the user has already alr-transformed the data.
#'
#' @param nus A 2 x 1 numeric vector with the smoothness coefficients for the Mat\'{e}rn covariance function, c(nu1, nu2).
#'
#'
#' @return \code{compositional_biwm_krig} returns a matrix with alr-transformed data for given krigging locations.
#'
#' @export
#'
#' @examples
#' set.seed(3)
#' n <- 40
#' coords <- matrix(runif(2*n), ncol = 2)
#' temp <- rnorm(2*n)
#' nus_vec <- c(0.5, 0.5)
#' nug_vec <- c(0,0)
#' S <- sigma_assembler_biwm(sigmas = c(1, 1), a = 2, rho = 0.5,
#' nus = nus_vec, coords_matrix = coords, combined = TRUE,nug_vec = nug_vec)
#' log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)
#'
#' informed_test <- fit_biwm(log_cd, coords, c(1, 1, 2, .5), nus_vec,nug_vec = nug_vec)
#' biwm_fit <- fit_biwm(log_cd, coords, c(.5, .5, 4, .6), nus_vec,nug_vec = nug_vec)
#'
#'x_grid <- seq(0,1,length.out = 10)
#'krig_locations <- as.matrix(expand.grid(x_grid,x_grid))
#'
#'compositional_biwm_krig(biwm_fit,krig_locations,coords,log_cd,nus_vec,nug_vec)
#'

compositional_biwm_krig <- function(biwm_fit,krig_locations,fit_locations,obs_matrix,nus,nug_vec){

  if(class(krig_locations) != "matrix" | class(fit_locations) != "matrix"){
    krig_locations <- as.matrix(krig_locations)
    fit_locations <- as.matrix(fit_locations)
  }

  sigma_hat <-
    sigma_assembler_biwm(
      sigmas = biwm_fit$theta[1:2],
      a = biwm_fit$theta[3]
      ,rho = biwm_fit$theta[4],
      coords_matrix = fit_locations,
      nus = nus,
      nug_vec = nug_vec,
      combined = TRUE
    )


  sigma_hat_with_new_loc <-
    cov_with_new_loc(
      loc_new = krig_locations,
      loc_obs = fit_locations,
      sigmas = biwm_fit$theta[1:2],
      a = biwm_fit$theta[3],
      rho = biwm_fit$theta[4],
      nus = nus,
      nug_vec = nug_vec,
      combined = TRUE
    )

  z_prediction <-
    setNames(
      as.data.frame(
        cbind(
          t(
            apply(
              matrix(
                sigma_hat_with_new_loc %*%
                  solve(sigma_hat,
                        as.vector(
                          t(
                            apply(
                              obs_matrix,
                              1,
                              function(x) x - biwm_fit$mu
                            )
                          )
                        )
                  ),
                ncol = 2
              ),
              1,
              function(x) x + biwm_fit$mu
            )
          ),
          krig_locations
        )
      ),
      c("process_01_pred","process_02_pred","coord_x","coord_y")
    )

  return(
    cbind(
      alr_inv(z_prediction[,c("process_01_pred","process_02_pred")]),
      z_prediction[,c("coord_x","coord_y")]
      )
  )
}
