#' Obtain Prediction Standard Error Estimate
#'
#'
#' @param krig_values A N x 2 numeric matrix of krig values not in the Simplex.
#'
#' @param biwm_fit An object of biwm class. Must be the result of a `fit_biwm` command.
#'
#' @param loc_new A N x 2 numeric matrix of coordinates for the unobserved data.
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
#' krig_vals <-
#'   compositional_biwm_krig(
#'     biwm_fit = generic_test,
#'     krig_locations = loc_new,
#'     fit_locations = coords,
#'     obs_matrix = log_cd,
#'     nus = nus,
#'     nug_vec = nug_vec,
#'     return_comp = FALSE
#'   )
#'
#' get_approx_alrInv_cond_covvar(krig_values = krig_vals,loc_new = loc_new, loc_obs = coords,biwm_fit = generic_test, nus = nus, nug_vec = nug_vec)
#'
#'
#'
#'

get_approx_alrInv_cond_covvar <- function(krig_values,loc_new,loc_obs,biwm_fit,nus,nug_vec){

  # ----- Aux Funs declaration -----

  extract_block_sigbar <- function(sig_bar){

    limits <- dim(sig_bar)/2

    block_sig11 <- sig_bar[1:limits[1],1:limits[2]]
    block_sig12 <- sig_bar[1:limits[1],(limits[2]+1):ncol(sig_bar)]
    block_sig21 <- sig_bar[(limits[1]+1):nrow(sig_bar),1:limits[2]]
    block_sig22 <- sig_bar[(limits[1]+1):nrow(sig_bar),(limits[2]+1):ncol(sig_bar)]

    return(
      list(
        block_sig11 = block_sig11,
        block_sig12 = block_sig12,
        block_sig21 = block_sig21,
        block_sig22 = block_sig22
      )
    )
  }



  extract_individual_sigbar <- function(iter,extract_block_sigbar_res){

    aux_mat <- matrix(0,2,2)

    aux_mat[1,1] <- extract_block_sigbar_res$block_sig11[iter,iter]
    aux_mat[1,2] <- extract_block_sigbar_res$block_sig12[iter,iter]
    aux_mat[2,1] <- extract_block_sigbar_res$block_sig21[iter,iter]
    aux_mat[2,2] <- extract_block_sigbar_res$block_sig22[iter,iter]

    return(
      aux_mat
    )
  }


  create_delta_grad_mat <- function(mu){

    #ffs
    if(is.data.frame(mu)){
      mu <- c(mu[[1]],mu[[2]])
    }

    dm_grad_mat <- matrix(0,2,3)

    common_denominator <- (1 + exp(mu[1]) + exp(mu[2]))^2

    dm_grad_mat[1,1] <- exp(mu[1]) * (1 + exp(mu[2]))
    dm_grad_mat[1,2] <- - exp(mu[1]) * exp(mu[2])
    dm_grad_mat[1,3] <- - exp(mu[1])

    dm_grad_mat[2,1] <- - exp(mu[1]) * exp(mu[2])
    dm_grad_mat[2,2] <- exp(mu[2]) * (1 + exp(mu[1]))
    dm_grad_mat[2,3] <- - exp(mu[2])



    return(
      dm_grad_mat/common_denominator
    )
  }



  # ------ Computation ------

  sigma_bar <- conditional_cov_mat(loc_new = loc_new,
                                      loc_obs = loc_obs,
                                      sigmas = biwm_fit$theta[1:2],
                                      a = biwm_fit$theta[3],
                                      rho = biwm_fit$theta[4],
                                      nus = nus,
                                      nug_vec = nug_vec)

  N <- nrow(sigma_bar) / 2

  sig_bar_blocks <- extract_block_sigbar(sig_bar = sigma_bar)

  approx_alrInv_cond_covvar <- matrix(0,N,3)

  for(i in 1:N){
    single_sig_bar <- extract_individual_sigbar(i,sig_bar_blocks)

    approx_mean <- create_delta_grad_mat(krig_values[i,])

    approx_alrInv_cond_covvar[i,] <-
      sqrt(diag(t(approx_mean) %*% single_sig_bar %*% approx_mean))
  }

  return(
    cbind(
      setNames(
        as.data.frame(approx_alrInv_cond_covvar),
        c("comp_01_sd","comp_02_sd","comp_03_sd")),
      setNames(
        as.data.frame(loc_new),
        c("coord_x","coord_y")
      )
    )
  )

}
