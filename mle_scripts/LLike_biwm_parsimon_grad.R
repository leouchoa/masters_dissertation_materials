source("../sigma_assembler_biwm.R")


####### -------  Assembling(organization purposes) function------------- 

matrices_assembler <- function(M_1,M_2,M_12){
  cbind(
    rbind(M_1,M_12),rbind(M_12,M_2)
  )
}

##### ----------- Matern derivative wrt aux function --------------

matern_deriv <- function(a,nu,coords_matrix){
  
  dist_vec <- as.vector(
    as.matrix(dist(coords_matrix))
  )
  # i'm following 'fields' package routine
  
  dist_vec[dist_vec == 0] <- 1e-10
  
  first_term <- ( 2^(nu - 1) * dist_vec^(nu) * a^(nu - 1) ) / 
    gamma(nu)
  
  snd_term <- 2 * nu * besselK(a * dist_vec,nu) -
    a * dist_vec * besselK(a * dist_vec, nu + 1)
  
  
  return(
    matrix(first_term*snd_term, 
           ncol = nrow(coords_matrix),
           nrow = nrow(coords_matrix)
    )
  )
  
}

####### ------  Log-likelihood general derivative form function ------- 

general_LLike_deriv <- function(cov_matrix,
                                grad_matrix, 
                                obs_vec,
                                ...){
  
  y <- solve(cov_matrix,obs_vec)
  
  -0.5 * (
    sum(diag(
      solve(cov_matrix, grad_matrix) 
    )) -
      as.double(t(y) %*% grad_matrix %*% y)
  )
  
}

####### ----------- main function : log-likelihood derivative ------

LLike_biwm_reduced_grad <- function(theta,
                                    nus,
                                    mu,
                                    coords_matrix,
                                    obs_vec,
                                    ...){
  
  other_params <- list(...)
  
  sigmas <- theta[1:2]
  as <- theta[3]
  rho <- theta[4]
  # mu <- theta[5]
  # nus <- theta[6:7]
  
  nus[3] <- mean(nus[1:2])
  if(length(as) == 1){
    as <- rep(as,3)
  }
  
  
  ##### ----- Non-gradient matrices -----
  
  #### Covariance matrix
  
  autocov_matrix <- sigma_assembler_biwm(sigmas = sigmas,
                                         as = as,
                                         rho = rho,
                                         nus = nus,
                                         coords_matrix = coords_matrix)
  
  #### Zero matrix
  
  zero_matrix <- matrix(0,
                        nrow = nrow(coords_matrix),
                        ncol = nrow(coords_matrix))
  
  ####### -------------  Gradient matrices ------------- 
  
  ##### Sigmas matrices
  
  sigmas_grad_matrices <- vector("list",2)
  
  sigmas_grad_matrices[[1]] <-  matrices_assembler(
    M_1 = matern_cov_wrapper(coords_matrix,
                             a = as[1],
                             nu = nus[1]),
    
    M_2 = zero_matrix,
    
    M_12 = (rho * sigmas[2] / sigmas[1] ) *
      matern_cov_wrapper(coords_matrix,
                         a = as[3],
                         nu = nus[3])
    )
  
  sigmas_grad_matrices[[2]] <-  matrices_assembler(
    
    M_1 = zero_matrix,
    
    M_2 = matern_cov_wrapper(coords_matrix,
                             a = as[2],
                             nu = nus[2]),
    
    M_12 = (rho * sigmas[1] / sigmas[2] ) *
      matern_cov_wrapper(coords_matrix,
                         a = as[3],
                         nu = nus[3])
  )
  
  ##### rho matrix
  
  rho_grad_matrix <- matrices_assembler(
    M_1 = zero_matrix,
    
    M_2 = zero_matrix,
    
    M_12 = sigmas[1] * sigmas[2] * matern_cov_wrapper(coords_matrix,
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
    
    M_12 = rho * sigmas[1] * sigmas[2] * matern_deriv(a = as[3], 
                                                      nu = nus[3],
                                                      coords_matrix = coords_matrix)
  )
  
  ###### ------ Eval of Gradients -----
  
  sigmas[1] <- general_LLike_deriv(cov_matrix = autocov_matrix,
                      grad_matrix = sigmas_grad_matrices[[1]],
                      obs_vec = obs_vec) 
  
  sigmas[2] <- general_LLike_deriv(cov_matrix = autocov_matrix,
                                   grad_matrix = sigmas_grad_matrices[[1]],
                                   obs_vec = obs_vec) 
  
  rho <- general_LLike_deriv(cov_matrix = autocov_matrix,
                            grad_matrix = rho_grad_matrix,
                            obs_vec = obs_vec)
  
  as <- general_LLike_deriv(cov_matrix = autocov_matrix,
                            grad_matrix = as_grad_matrix,
                            obs_vec = obs_vec)
  
  
  return(c(
    sigmas, a = as, rho = rho
  ))
}
