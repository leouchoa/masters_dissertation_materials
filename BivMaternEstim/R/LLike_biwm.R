#' Log-likelihood of the Bivariate (Wittle-)Matern Model.
#'
#' @param theta A 4 x 1 numeric vector, the order of parameters is: theta = c(sigma_1^2, sigma_2^2, a, rho).
#' @param nus A 2 x 1 numeric vector with the smoothness coefficients for the Mat\'{e}rn covariance function, c(nu1, nu2).
#' @param mu A 2 x 1 numeric vector of means c(mu1, mu2).
#' @param coords_matrix A n x 2 numeric matrix of coordinates for the data.
#' @param obs_matrix A n x 2 numeric matrix of data. We assume that the user has already log-transformed the data.
#'
#' @return log-likelihood
#'
#' @export
LLike_biwm <- function(theta,
                       nus,
                       mu,
                       coords_matrix,
                       obs_matrix){

  sigmas <- theta[1:2]
  as <- theta[3]
  rho <- theta[4]

  if(length(as) == 1){
    as <- rep(as, 3)
  }

  autocov_matrix <- sigma_assembler_biwm(sigmas = sigmas,
                                         a = as,
                                         rho = rho,
                                         nus = nus,
                                         coords_matrix = coords_matrix)

  autocov_matrix_inv <- block_inv(autocov_matrix)

  Z <- t(apply(obs_matrix, 1, function(x) x - mu))

  as.double(
    -0.5 * (
      block_log_det(autocov_matrix) +
        t(Z[,1]) %*% autocov_matrix_inv[[1]] %*% as.matrix(Z[,1]) +
        2 * t(Z[,1]) %*% autocov_matrix_inv[[2]] %*% as.matrix(Z[,1]) +
        t(Z[,2]) %*% autocov_matrix_inv[[3]] %*% as.matrix(Z[,1]) +
        length(obs_matrix) * log(2*pi)
    ))


  # as.double(
  #   -1/2 *
  #     (
  #       log(det(autocov_matrix)) +
  #         crossprod(Z, solve(autocov_matrix, Z)) +
  #         length(obs_matrix) * log(2*pi)
  #     )
  # )

}


# ---- Sanity Check ----
#
# aux_assembler_inv <- function(matrix_list){
#   cbind(
#     rbind(matrix_list[[1]], t(matrix_list[[2]])),
#           rbind(matrix_list[[2]], matrix_list[[3]])
#     )
# }
#
# aux_assembler_autocov <- function(matrix_list_2){
#   cbind(
#     rbind(matrix_list_2[[1]], matrix_list_2[[3]]),
#     rbind(matrix_list_2[[3]], matrix_list_2[[2]])
#   )
# }
#
# autocov_matrix_old <- aux_assembler_autocov(autocov_matrix)
#
# all.equal(log(det(autocov_matrix_old)),
#           block_log_det(autocov_matrix))
#
# S_inv_old <- solve(autocov_matrix_old)
# S_inv_new <- aux_assembler_inv(autocov_matrix_inv)
#
# all.equal(S_inv_new,S_inv_old) #hmmm
#
#
# kernel_new <- t(Z[,1]) %*% autocov_matrix_inv[[1]] %*% Z[,1] +
#   2 * t(Z[,1]) %*% autocov_matrix_inv[[2]] %*% Z[,2] +
#   t(Z[,2]) %*% autocov_matrix_inv[[3]] %*% Z[,2]
#
# kernel_old <- crossprod(as.numeric(Z), solve(autocov_matrix_old, as.numeric(Z)))
#
# kernel_test_new <- t(as.numeric(Z)) %*% S_inv_new %*% as.numeric(Z)
#
# log_lik_new_way <- -0.5 * (
#   block_log_det(autocov_matrix) +
#     t(Z[,1]) %*% autocov_matrix_inv[[1]] %*% Z[,1] +
#     2 * t(Z[,1]) %*% autocov_matrix_inv[[2]] %*% Z[,2] +
#     t(Z[,2]) %*% autocov_matrix_inv[[3]] %*% Z[,2] +
#     length(obs_matrix) * log(2*pi)
#
# )
#
# log_lik_old_way <-
#   as.double(
#     -1/2 *
#       (
#         log(det(autocov_matrix_old)) +
#           crossprod(as.numeric(Z), solve(autocov_matrix_old, as.numeric(Z))) +
#           length(obs_matrix) * log(2*pi)
#       )
#   )
#
# log_lik_new_inv_old_mahala <-
#
#   as.double(
#     -1/2 *
#       (
#         log(det(autocov_matrix_old)) +
#           t(as.numeric(Z)) %*% S_inv_new %*% as.numeric(Z)
#         +
#           length(obs_matrix) * log(2*pi)
#       )
#
#   )
