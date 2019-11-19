# equivalence between RandomFields notation and their article.... ffs
# 1. r = h 
# 2. s_ij = 1/a_ij
# 3. c_ii = sigma_ii

source("../LLike_biwm.R")
source("../LLike_biwm_parsimon_grad.R")

# rho_test <- function(rho,nus){
#   if(abs(rho) > ( sqrt(nus[1]*nus[2]) )/mean(nus)) 
#     print("Not a valid covariance")
#   else print("Valid covariance")
# }

library(RandomFields)
x_test <- y_test <- seq(-10, 10, length.out = 10)
model <- RMbiwm(
  nudiag = c(1, 2), 
  nured = 1.5, 
  rhored = 0.6, 
  cdiag = c(3, 6),
  s = c(2, 2, 2))

simu_parsimon <- RFsimulate(model, x_test, y_test)
# X11();plot(simu_parsimon)


true_params <- c(sigmas = c(3,6),
                 as = 0.5,
                 rho = 0.6
                 ) 

init_params = c(sigmas = c(3,6),
                as = 0.5,
                rho = 0.6
                )

# init_params_test <- c(3,6,0.5,0.6)

# result_parsimon <- optim(par = init_params,
#                 control = list(fnscale = -1),
#                 fn = LLike_biwm,
#                 method = "L-BFGS-B",
#                 lower = c(0.001,0.001,-1,0.001),
#                 upper = c(Inf,Inf,1,Inf),
#                 nus = c(1,1.5),
#                 parsimonious = TRUE,
#                 coords_matrix = coordinates(simu_parsimon),
#                 obs_vec = c(simu_parsimon$variable1,simu_parsimon$variable2)
#                 )

# saveRDS(result,"optim_result.rds")
# readRDS(result,"optim_result.rds")



result_parsimon_wgrad <- optim(par = init_params,
                               control = list(fnscale = -1,
                                              trace = 4),
                               fn = LLike_biwm,
                               gr = LLike_biwm_reduced_grad,
                               method = "L-BFGS-B",
                               lower = c(0.001,0.001,0.001,-1),
                               upper = c(Inf,Inf,Inf,1),
                               nus = c(1,2),
                               mu = 0,
                               coords_matrix = coordinates(simu_parsimon),
                               obs_vec = c(simu_parsimon$variable1,simu_parsimon$variable2)
)
 
 # Error in optim(par = init_params, control = list(fnscale = -1), fn = LLike_biwm,  : 
 #                  L-BFGS-B needs finite values of 'fn'
 # 
 # wtf???? The code below doesn't return NA..... 

 LLike_biwm(init_params,
            c(1,2),
            0,
            coordinates(simu_parsimon),
            c(simu_parsimon$variable1,simu_parsimon$variable2))

 LLike_biwm_reduced_grad(init_params,
            c(1,2),
            0,
            coordinates(simu_parsimon),
            c(simu_parsimon$variable1,simu_parsimon$variable2))
