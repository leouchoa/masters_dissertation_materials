source("../LLike_biwm_parsimon_grad.R")

library(RandomFields)
x_test <- y_test <- seq(-10, 10, length.out = 10)
model <- RMbiwm(
  nudiag = c(1, 2), 
  nured = 1.5, 
  rhored = 0.6, 
  cdiag = c(3, 6),
  s = c(2, 2, 2))

simu_parsimon <- RFsimulate(model, x_test, y_test)

####### ---- Testing of construction matrices ------

# true_params <- c(sigmas = c(3,6),
#                  rho = 0.6,
#                  as = 0.5,
#                  nus = c(1,1.5)
# ) 

init_params <- c(sigmas_1 = 3,
                    sigmas_2 = 6,
                    as = 0.5,
                    rho = 0.6,
                    mu = 1, #not true
                    nus_1 = 1,
                    nus_2 = 1.5
) 

coords_matrix = coordinates(simu_parsimon)
obs_vec = c(simu_parsimon$variable1,simu_parsimon$variable2)

theta <- init_params
sigmas <- theta[1:2]
as <- theta[3]
rho <- theta[4]
mu <- theta[5]
nus <- theta[6:7]
nus[3] <- mean(nus[1:2])
if(length(as) == 1){
  as <- rep(as,3)
}

  ##### ----- Non-gradient matrices -----
  
  #### Covariance matrix
  
  autocov_matrix <- sigma_assembler_biwm(coords_matrix = coords_matrix,
                                         sigmas = sigmas,
                                         as = as,
                                         rho = rho,
                                         nus = nus)
  
  # image(autocov_matrix)
  # ok passed
  
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
  
  # image(sigmas_grad_matrices[[1]])
  # image(sigmas_grad_matrices[[2]]) # strange af, wth
  # ok passed
  
  ##### rho matrix
  
  rho_grad_matrix <- matrices_assembler(
    M_1 = zero_matrix,
    
    M_2 = zero_matrix,
    
    M_12 = sigmas[1] * sigmas[2] * matern_cov_wrapper(coords_matrix,
                                                      a = as[1],
                                                      nu = nus[1])
    
  )
  # image(rho_grad_matrix)
  # ok passed
  
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
  
  # image(as_grad_matrix) #WTFFFFFFFFFFFFFFFFFFFFFFF
  # fk
  
  ###### ------ Eval of Gradients -----
  
  sigmas[1] <- general_LLike_deriv(cov_matrix = autocov_matrix,
                                   grad_matrix = sigmas_grad_matrices[[1]],
                                   obs_vec = obs_vec) 
  
  sigmas[2] <- general_LLike_deriv(cov_matrix = autocov_matrix,
                                   grad_matrix = sigmas_grad_matrices[[2]],
                                   obs_vec = obs_vec) 
  
  # sigmas # holy fk, fk
  
  rho <- general_LLike_deriv(cov_matrix = autocov_matrix,
                             grad_matrix = rho_grad_matrix,
                             obs_vec = obs_vec)
  # ok passed
  
  as <- general_LLike_deriv(cov_matrix = autocov_matrix,
                            grad_matrix = as_grad_matrix,
                            obs_vec = obs_vec)
  
  # ---fk (again)--
  # update: got to solve this by following fields package in 0 zero distance
  # besselK 
  
  # y <- solve(autocov_matrix,obs_vec)
  # 
  # sum(diag(
  #   solve(autocov_matrix) %*% as_grad_matrix
  # )) -
  #   as.double(t(y) %*% as_grad_matrix %*% y)
  
  
  c(sigmas[1], sigmas[2], as, rho)

  
  sigmas = c(3,6)
  rho = 0.6
  as = rep(0.5,3)
  nus = c(1,1.5);nus[3] <- mean(nus)
  
  coords_matrix = coordinates(simu_parsimon)
  obs_vec = c(simu_parsimon$variable1,simu_parsimon$variable2)
    
  res <- LLike_biwm_reduced_grad(sigmas = sigmas,
                          rho = rho,
                          as = as,
                          nus = nus,
                          coords_matrix = coords_matrix,
                          obs_vec = obs_vec)
  
  # both agreed!

