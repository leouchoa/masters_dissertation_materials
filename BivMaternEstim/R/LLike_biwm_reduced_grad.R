LLike_biwm_reduced_grad <- function(theta,
                                    nus,
                                    mu,
                                    coords_matrix,
                                    obs_vec){

  Z <- as.numeric(t(apply(obs_vec, 1, function(x) x - mu)))

  sigmas <- theta[1:2]
  as <- theta[3]
  rho <- theta[4]

  nus[3] <- mean(nus[1:2])

  if(length(as) == 1){
    as <- rep(as,3)
  }

  #### Covariance matrix

  autocov_matrix <- sigma_assembler_biwm(sigmas = sigmas,
                                         a = as,
                                         rho = rho,
                                         nus = nus,
                                         coords_matrix = coords_matrix)

  #### Zero matrix

  zero_matrix <- matrix(0,
                        nrow = nrow(coords_matrix),
                        ncol = nrow(coords_matrix))

  ####### -------------  Gradient matrices -------------

  d <- dist(coords_matrix)

  ##### Sigmas matrices

  sigmas_grad_matrices <- vector("list",2)

  sigmas_grad_matrices[[1]] <-  matrices_assembler(
    M_1 = matern_cov_wrapper(d,
                             a = as[1],
                             nu = nus[1]),

    M_2 = zero_matrix,

    M_12 = 0.5*(rho * sqrt(sigmas[2] / sigmas[1]) ) *
      matern_cov_wrapper(d,
                         a = as[3],
                         nu = nus[3])
  )

  sigmas_grad_matrices[[2]] <-  matrices_assembler(

    M_1 = zero_matrix,

    M_2 = matern_cov_wrapper(d,
                             a = as[2],
                             nu = nus[2]),

    M_12 = 0.5*(rho * sqrt(sigmas[1] / sigmas[2]) ) *
      matern_cov_wrapper(d,
                         a = as[3],
                         nu = nus[3])
  )

  ##### rho matrix

  rho_grad_matrix <- matrices_assembler(
    M_1 = zero_matrix,

    M_2 = zero_matrix,

    M_12 = sqrt(sigmas[1] * sigmas[2]) * matern_cov_wrapper(d,
                                                            a = as[1],
                                                            nu = nus[1])

  )

  ##### As matrix


  as_grad_matrix <- matrices_assembler(

    M_1 = sigmas[1] * matern_deriv(a = as[1],
                                   nu = nus[1],
                                   coords_matrix = coords_matrix),

    M_2 = sigmas[2] * matern_deriv(a = as[2],
                                   nu = nus[2],
                                   coords_matrix = coords_matrix),

    M_12 = rho * sqrt(sigmas[1] * sigmas[2]) * matern_deriv(a = as[3],
                                                            nu = nus[3],
                                                            coords_matrix = coords_matrix)
  )

  ###### ------ Eval of Gradients -----

  sigmas[1] <- general_LLike_deriv(cov_matrix = autocov_matrix,
                                   grad_matrix = sigmas_grad_matrices[[1]],
                                   obs_vec = Z)

  sigmas[2] <- general_LLike_deriv(cov_matrix = autocov_matrix,
                                   grad_matrix = sigmas_grad_matrices[[1]],
                                   obs_vec = Z)

  rho <- general_LLike_deriv(cov_matrix = autocov_matrix,
                             grad_matrix = rho_grad_matrix,
                             obs_vec = Z)

  as <- general_LLike_deriv(cov_matrix = autocov_matrix,
                            grad_matrix = as_grad_matrix,
                            obs_vec = Z)


  return(c(sigmas, a = as, rho = rho))
}
