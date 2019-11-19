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

n_coord = coordinates(simu_parsimon)
obs_vec = c(simu_parsimon$variable1,simu_parsimon$variable2)


init_params <- c(sigmas_1 = 3,
           sigmas_2 = 6,
           as = 0.5,
           rho = 0.6)

mu = 1 #not true
nus = c(1,2)


LLike_biwm_reduced_grad(init_params,
                        nus = nus,
                        mu = 1,
                        coords_matrix = n_coord,
                        obs_vec = obs_vec)

