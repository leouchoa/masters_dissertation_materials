block_LLike_biwm_grad <- function(theta,
                                    nus,
                                    mu,
                                    coords_matrix,
                                    obs_matrix){



  Z <- as.numeric(t(apply(obs_matrix, 1, function(x) x - mu)))

  sigmas <- theta[1:2]
  as <- theta[3]
  rho <- theta[4]

  nus[3] <- mean(nus[1:2])

  if(length(as) == 1){
    as <- rep(as,3)
  }

  # ---- Covariance matrix and its inverse ----

  autocov_matrix <- sigma_assembler_biwm(sigmas = sigmas,
                                         a = as,
                                         rho = rho,
                                         nus = nus,
                                         coords_matrix = coords_matrix)

  autocov_matrix_inv <- block_inv(autocov_matrix)

  y_1 <- autocov_matrix_inv$C_11_star %*% obs_matrix[,1] +
               autocov_matrix_inv$C_12_star %*% obs_matrix[,2]

  y_2 <- t(autocov_matrix_inv$C_12_star) %*% obs_matrix[,1] +
               autocov_matrix_inv$C_22_star %*% obs_matrix[,2]

  # ----- Grad of sigma^2_1 -----
  sigma_1_tr <-
    (
      tr(autocov_matrix_inv$C_11_star %*% autocov_matrix$C_11) +
        tr(autocov_matrix_inv$C_12_star %*% autocov_matrix$C_12)
    ) / sigmas[1]

  sigma_1_qf <-
    as.double(t(y_1) %*%
                (autocov_matrix$C_11 %*% y_1 +
                   autocov_matrix$C_12 %*% y_2))/ sigmas[1]

  sigma_1_grad <- -0.5*(sigma_1_tr - sigma_1_qf)


  # ----- Grad of sigma^2_2 -----


  sigma_2_tr <-
    (
      tr(autocov_matrix_inv$C_22_star %*% (autocov_matrix$C_22)) +
        tr(autocov_matrix_inv$C_12_star %*% (autocov_matrix$C_12))
    ) / sigmas[2]

  sigma_2_qf <-
    as.double(
      (t(y_2) %*% autocov_matrix$C_22 +
         t(y_1) %*% autocov_matrix$C_12) %*% y_2) / sigmas[2]

  sigma_2_grad <- -0.5*(sigma_2_tr - sigma_2_qf)


  # ----- Grad of rho -----

  M_nu_3 <- matern_cov_wrapper(coords_dist = dist(coords_matrix),
                               a = as[3],
                               nu = nus[3])

  rho_tr <- 2*sqrt(sigmas[1]*sigmas[2]) *
    tr(
      autocov_matrix_inv$C_12_star %*% M_nu_3
    )

  rho_qf <- as.double(
    2*sqrt(sigmas[1]*sigmas[2]) *
      t(y_1) %*% M_nu_3 %*% y_2)

  rho_grad <- -0.5*(rho_tr - rho_qf)

  # ----- Grad of a -----

  M_dash_nu_1 <- matern_deriv(a = as[1], nu = nus[1], coords_matrix = coords_matrix)

  M_dash_nu_2 <- matern_deriv(a = as[2], nu = nus[2], coords_matrix = coords_matrix)

  M_dash_nu_3 <- matern_deriv(a = as[3], nu = nus[3], coords_matrix = coords_matrix)

  a_tr <- sigmas[1] * tr(autocov_matrix_inv$C_11_star %*% M_dash_nu_1) +
    sigmas[2] * tr(autocov_matrix_inv$C_22_star %*% M_dash_nu_2) +
    2*rho*sqrt(sigmas[1]*sigmas[2]) * tr(autocov_matrix_inv$C_12_star %*% M_dash_nu_3)

  a_qf <- as.double(
    sigmas[1] * t(y_1) %*% M_dash_nu_1 %*% y_1 +
      sigmas[2] * t(y_2) %*% M_dash_nu_2 %*% y_2 +
      2*rho*sqrt(sigmas[1]*sigmas[2])  * t(y_1) %*% M_dash_nu_3 %*% y_2)

  a_grad <- -0.5*(a_tr - a_qf)


    return(c(
      sigmas = c(sigma_1_grad,sigma_2_grad),
      a = a_grad,
      rho = rho_grad
      ))

}


# # ---- Usage ----
#
# source("block_inv.R")
# source("utils.R")
# source("matern_deriv.R")
# source("matern_cov_wrapper.R")
# source("sigma_assembler_biwm.R")
#
# set.seed(1233)
# n <- 40
# theta_test <- c(1,1,2,0.5)
# coords_test <- matrix(runif(2*n), ncol = 2)
# nus_test <- rep(0.5,2)
# S <- sigma_assembler_biwm(sigmas = c(1, 1), a = 2, rho = 0.5, nus = nus_test, coords_matrix = coords_test,combined = TRUE)
# temp <- rnorm(2*n)
# log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)
#
# new_grad <- block_LLike_biwm_grad(theta_test,nus_test,colMeans(log_cd),coords_test,log_cd)
#
# source("../old_code/matrices_assembler.R")
# source("../old_code/general_LLike_deriv.R")
# source("../old_code/LLike_biwm_reduced_grad.R")
# source("../old_code/sigma_assembler_biwm.R")
#
# old_grad <- LLike_biwm_reduced_grad(theta_test,nus_test,colMeans(log_cd),coords_test,log_cd)
#
