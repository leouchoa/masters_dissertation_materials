get_fisher_info <- function(theta_vec,nus,coords_matrix){
  
  tr <- function(x){sum(diag(x))}
  
  # ------ Params Extraction -------  
  
  sigmas <- theta_vec[1:2]
  as <- rep(theta_vec[3],3)
  rho <- theta_vec[4]
  nus[3] <- mean(nus)
  
  d <- dist(coords_matrix)
  
  
  # ------ Inverse Calculation -------  
  
  sigma_matrix <- BivMaternEstim:::sigma_assembler_biwm(sigmas,
                                                        as[1],
                                                        rho,
                                                        nus,
                                                        coords_matrix,
                                                        FALSE
  )
  
  sigma_matrix_inv <- BivMaternEstim:::block_inv(sigma_matrix)
  C_11_st <- sigma_matrix_inv$C_11_star
  C_12_st <- sigma_matrix_inv$C_12_star
  C_22_st <- sigma_matrix_inv$C_22_star
  # -------------- Derivative matrices ---------------
  
  
  #' # Pregame
  #' 
  #' 
  #' ## Creating Matern Matrices
  #' 
  #' 
  #' 
  
  
  M_1 <- BivMaternEstim:::matern_cov_wrapper(d,a = as[1],nu = nus[1])
  M_2 <- BivMaternEstim:::matern_cov_wrapper(d,a = as[2],nu = nus[2])
  M_3 <- BivMaternEstim:::matern_cov_wrapper(d,a = as[3],nu = nus[3])
  
  #' ## Creating Matern Derivative Matrices
  #' 
  #' 
  
  M_dash_nu_1 <- BivMaternEstim:::matern_deriv(a = as[1], 
                                               nu = nus[1], 
                                               coords_matrix = coords_matrix)
  
  M_dash_nu_2 <- BivMaternEstim:::matern_deriv(a = as[2], 
                                               nu = nus[2], 
                                               coords_matrix = coords_matrix)
  
  M_dash_nu_3 <- BivMaternEstim:::matern_deriv(a = as[3], 
                                               nu = nus[3], 
                                               coords_matrix = coords_matrix)
  
  #' # Caso em que $\theta = \sigma^2_1$ 
  #' 
  #' $B = \Sigma^{-1} \frac{\partial \Sigma}{\partial \sigma^2_1}$
  #' 
  #' 
  #' 
  #' 
  
  B_11_sigma1 = C_11_st %*% M_1 +   rho * sqrt( sigmas[2]/(2*sigmas[1]) )* C_12_st %*% M_3
  B_21_sigma1 = t(C_12_st) %*% M_1 +  rho * sqrt( sigmas[2]/(2*sigmas[1]) ) * C_22_st %*% M_3
  B_12_sigma1 = rho * sqrt( sigmas[2]/(2*sigmas[1]) ) * C_11_st %*% M_3
  B_22_sigma1 = rho * sqrt( sigmas[2]/(2*sigmas[1]) ) * t(C_12_st) %*% M_3
  
  
  #' # Caso em que $\theta = \sigma^2_1$ 
  #' 
  #' $B = \Sigma^{-1} \frac{\partial \Sigma}{\partial \sigma^2_1}$
  #' 
  #' 
  #' 
  #' 
  
  B_11_sigma2 = rho * sqrt( sigmas[1]/(2*sigmas[2]) ) * C_12_st %*% M_3
  B_21_sigma2 = rho * sqrt( sigmas[1]/(2*sigmas[2]) ) * C_22_st %*% M_3
  B_12_sigma2 = rho * sqrt( sigmas[1]/(2*sigmas[2]) ) * C_11_st %*% M_3 + C_12_st %*% M_2
  B_22_sigma2 = rho * sqrt( sigmas[1]/(2*sigmas[2]) ) * t(C_12_st) %*% M_3 + C_22_st %*% M_2
  
  
  #' # Caso em que $\theta = \rho$ 
  #' 
  #' $B = \Sigma^{-1} \frac{\partial \Sigma}{\partial \rho}$
  #' 
  #' 
  #' 
  #' 
  
  B_11_rho = sqrt( sigmas[1]* sigmas[2] ) * C_12_st %*% M_3
  B_21_rho = sqrt( sigmas[1]* sigmas[2] ) * C_22_st %*% M_3
  B_12_rho = sqrt( sigmas[1]* sigmas[2] ) * C_11_st %*% M_3
  B_22_rho = sqrt( sigmas[1]* sigmas[2] ) * t(C_12_st) %*% M_3
  
  
  
  #' # Caso em que $\theta = a$ 
  #' 
  #' $B = \Sigma^{-1} \frac{\partial \Sigma}{\partial \rho}$
  #' 
  #' 
  #' 
  #' 
  
  B_11_a = sigmas[1] * C_11_st %*% M_dash_nu_1 + 
    rho *sqrt( sigmas[1]* sigmas[2] ) * C_12_st %*% M_dash_nu_3
  
  B_21_a = sigmas[1] * t(C_12_st) %*% M_dash_nu_1 + 
    rho *sqrt( sigmas[1]* sigmas[2] ) * C_22_st %*% M_dash_nu_3
  
  B_12_a = sigmas[2] * C_12_st %*% M_dash_nu_2 + 
    rho *sqrt( sigmas[1]* sigmas[2] ) * C_11_st %*% M_dash_nu_3
  
  B_22_a = sigmas[2] * C_22_st %*% M_dash_nu_2 + 
    rho *sqrt( sigmas[1]* sigmas[2] ) * t(C_12_st) %*% M_dash_nu_3
  
  
  # -------------- Fisher Information ---------------
  
  #' # Fisher information first line
  #' 
  #' 
  #' 
  
  #' ## Fisher information for $\sigma^2_1$
  
  F_11 <- 0.25 *(
    2 * tr(
      B_11_sigma1 %*% B_11_sigma1  + B_12_sigma1 %*% B_21_sigma1
    ) + 
      2 * tr(
        B_21_sigma1 %*% B_12_sigma1  + B_22_sigma1 %*% B_22_sigma1
      ) - 
      (tr( B_11_sigma1 ) + tr( B_22_sigma1 ))^2
  )
  
  #' ## Fisher information for $\sigma^2_1$ and $\sigma^2_2$
  
  F_12 <- 0.25 *(
    2 * tr(
      B_11_sigma1 %*% B_11_sigma2  + B_12_sigma1 %*% B_21_sigma2
    ) + 
      2 * tr(
        B_21_sigma1 %*% B_12_sigma2  + B_22_sigma1 %*% B_22_sigma2
      ) - 
      (tr( B_11_sigma1 ) + tr( B_22_sigma1 )) * 
      (tr( B_11_sigma2 ) + tr( B_22_sigma2 ))
  )
  
  #' ## Fisher information for $\sigma^2_1$ and $\rho$
  
  
  F_13 <- 0.25 *(
    2 * tr(
      B_11_sigma1 %*% B_11_rho  + B_12_sigma1 %*% B_21_rho
    ) + 
      2 * tr(
        B_21_sigma1 %*% B_12_rho  + B_22_sigma1 %*% B_22_rho
      ) - 
      (tr( B_11_sigma1 ) + tr( B_22_sigma1 )) * 
      (tr( B_11_rho ) + tr( B_22_rho ))
  )  
  
  
  #' ## Fisher information for $\sigma^2_1$ and $a$
  
  
  F_14 <- 0.25 *(
    2 * tr(
      B_11_sigma1 %*% B_11_a  + B_12_sigma1 %*% B_21_a
    ) + 
      2 * tr(
        B_21_sigma1 %*% B_12_a  + B_22_sigma1 %*% B_22_a
      ) - 
      (tr( B_11_sigma1 ) + tr( B_22_sigma1 )) * 
      (tr( B_11_a ) + tr( B_22_a ))
  )  
  
  #' # Fisher information second line
  #' 
  #' 
  #' ## Fisher information for $\sigma^2_2$
  
  F_22 <- 0.25 *(
    2 * tr(
      B_11_sigma2 %*% B_11_sigma2  + B_12_sigma2 %*% B_21_sigma2
    ) + 
      2 * tr(
        B_21_sigma2 %*% B_12_sigma2  + B_22_sigma2 %*% B_22_sigma2
      ) - 
      (tr( B_11_sigma2 ) + tr( B_22_sigma2 ))^2
  )
  
  
  #' ## Fisher information for $\sigma^2_2$ and $\rho$
  #' 
  
  
  F_23 <- 0.25 *(
    2 * tr(
      B_11_sigma2 %*% B_11_rho  + B_12_sigma2 %*% B_21_rho
    ) + 
      2 * tr(
        B_21_sigma2 %*% B_12_rho  + B_22_sigma2 %*% B_22_rho
      ) - 
      (tr( B_11_sigma2 ) + tr( B_22_sigma2 )) * 
      (tr( B_11_rho ) + tr( B_22_rho ))
  ) 
  
  F_24 <- 0.25 *(
    2 * tr(
      B_11_sigma2 %*% B_11_a  + B_12_sigma2 %*% B_21_a
    ) + 
      2 * tr(
        B_21_sigma2 %*% B_12_a  + B_22_sigma2 %*% B_22_a
      ) - 
      (tr( B_11_sigma2 ) + tr( B_22_sigma2 )) * 
      (tr( B_11_a ) + tr( B_22_a ))
  ) 
  
  #' # Fisher information third line
  #' 
  #' 
  #' 
  #' ## Fisher information for $\rho$
  
  F_33 <- 0.25 *(
    2 * tr(
      B_11_rho %*% B_11_rho  + B_12_rho %*% B_21_rho
    ) + 
      2 * tr(
        B_21_rho %*% B_12_rho  + B_22_rho %*% B_22_rho
      ) - 
      (tr( B_11_rho ) + tr( B_22_rho ))^2
  )
  
  #' 
  #' ## Fisher information for $\rho$ and $a$
  #' 
  #' 
  #' 
  
  F_34 <- 0.25 *(
    2 * tr(
      B_11_rho %*% B_11_a  + B_12_rho %*% B_21_a
    ) + 
      2 * tr(
        B_21_rho %*% B_12_a  + B_22_rho %*% B_22_a
      ) - 
      (tr( B_11_rho ) + tr( B_22_rho )) * 
      (tr( B_11_a ) + tr( B_22_a ))
  ) 
  
  
  #' # Fisher information fourth line
  #' 
  #' 
  
  F_44 <- 0.25 *(
    2 * tr(
      B_11_a %*% B_11_a  + B_12_a %*% B_21_a
    ) + 
      2 * tr(
        B_21_a %*% B_12_a  + B_22_a %*% B_22_a
      ) - 
      (tr( B_11_a ) + tr( B_22_a ))^2
  )
  
  # --------------- Assembling the expected values -----------
  
  fisher_info <- matrix(0,4,4)
  
  fisher_info[1,1] <- F_11
  fisher_info[1,2] <- F_12
  fisher_info[1,3] <- F_13
  fisher_info[1,4] <- F_14
  
  fisher_info[2,1] <- F_12
  fisher_info[2,2] <- F_22
  fisher_info[2,3] <- F_23
  fisher_info[2,4] <- F_24
  
  fisher_info[3,1] <- F_13
  fisher_info[3,2] <- F_23
  fisher_info[3,3] <- F_33
  fisher_info[3,4] <- F_34
  
  fisher_info[4,1] <- F_14
  fisher_info[4,2] <- F_24
  fisher_info[4,3] <- F_34
  fisher_info[4,4] <- F_44
  
  fisher_info
}
