source("R/matern_cov_wrapper.R")

sigma_assembler_biwm <- function(sigmas,
                                 as,
                                 rho,
                                 nus,
                                 coords_matrix
                                 ){

  M_1 <- sigmas[1]* matern_cov_wrapper(coords_matrix,
                                       a = as[1],
                                       nu = nus[1])

  M_2 <- sigmas[2]* matern_cov_wrapper(coords_matrix,
                                       a = as[2],
                                       nu = nus[2])

  M_12 <- rho * sigmas[1] * sigmas[2] * matern_cov_wrapper(coords_matrix,
                                                         a = as[3],
                                                         nu = nus[3])
  # return(list(M_1 = M1,
  #             M_2 = M2,
  #             M_12 = M12))

  return(
    cbind(
      rbind(M_1,M_12),rbind(M_12,M_2)
    )
  )
}
