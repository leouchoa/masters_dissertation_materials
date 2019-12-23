sigma_assembler_biwm <- function(sigmas, as, rho, nus, coords_matrix){

  d <- dist(coords_matrix)

  M_1 <- sigmas[1] * matern_cov_wrapper(d,
                                        a = as[1],
                                        nu = nus[1])

  M_2 <- sigmas[2] * matern_cov_wrapper(d,
                                        a = as[2],
                                        nu = nus[2])

  M_12 <- rho * sqrt(sigmas[1] * sigmas[2]) * matern_cov_wrapper(d,
                                                           a = as[3],
                                                           nu = nus[3])

  ret <- cbind(rbind(M_1, t(M_12)), rbind(M_12, M_2))

  return(ret)
}
