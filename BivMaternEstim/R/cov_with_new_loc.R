cov_with_new_loc <- function(loc_new,loc_obs,sigmas, a, rho, nus, nug_vec, combined = FALSE){


  # NEED TO DO POINTWISE DISTANCE
  d <- as.matrix(dist(rbind(loc_obs,loc_new)))
  #
  # # I don't know if next line is necessary. The point is that the krigging for the observed locations is exact, although the standard error map is not
  #
  d <- d[ (nrow(loc_obs) + 1) : nrow(d),  1:nrow(loc_obs) ]

  # ---- LEFT OFF HERE -----

  if(length(nus) == 2) nus[3] <- mean(nus[1:2])

  M_1 <- sigmas[1] * matern_cov_wrapper_krig(d,
                                        a = a,
                                        nu = nus[1])
  diag(M_1) <- diag(M_1)

  M_2 <- sigmas[2] * matern_cov_wrapper_krig(d,
                                        a = a,
                                        nu = nus[2])

  diag(M_2) <- diag(M_2)

  M_12 <- rho * sqrt(sigmas[1] * sigmas[2]) * matern_cov_wrapper_krig(d,
                                                                 a = a,
                                                                 nu = nus[3])


  return(
    cbind(rbind(M_1, M_12), rbind(M_12, M_2))
    )

}
