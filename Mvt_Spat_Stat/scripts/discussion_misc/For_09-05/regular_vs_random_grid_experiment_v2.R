# ---- 01 global settings ----
library(BivMaternEstim)

set.seed(1234)
n_points <- 7
temp <- rnorm(2* (n_points^2))
nus_vec <- c(0.5, 0.5)

true_param <- c(1.5,0.7,1.3,0.4)
initial_guess_1 <-  c(0.5,1.2,4,-0.2)
initial_guess_2 <-  c(0.5,1.2,0.5,0.2)
initial_guess_3 <-  c(0.5,1.2,0.1,0.8)


# ---- 02 regular grid ----


x <- seq(-1,1,length.out = n_points)
coords_regular <- expand.grid(x,x)


S_regular <- sigma_assembler_biwm(sigmas = true_param[1:2], a = true_param[3], rho = true_param[4], nus = nus_vec, coords_matrix = coords_regular, combined = TRUE)

log_cd_regular <- matrix(temp%*%chol(S_regular) + rep(c(1,2), each = n_points), ncol = 2)

informed_test_regular <- fit_biwm(log_cd_regular, coords_regular, true_param, nus_vec)

generic_test_regular <- fit_biwm(log_cd_regular, coords_regular, initial_guess_1, nus_vec)



exp_res_regular <- 
  as.data.frame(
    rbind(
      true_param = round(true_param,5),
      initial_guess = round(initial_guess_1,5),
      estimated_generic = round(generic_test_regular$theta,5),
      estimated_informed = round(informed_test_regular$theta,5)
    )
  )




# ---- 03  random grid ----


coords_random <- matrix(runif(2 * (n_points^2) ,-1,1), ncol = 2)

S_random <- sigma_assembler_biwm(sigmas = true_param[1:2], a = true_param[3], rho = true_param[4], nus = nus_vec, coords_matrix = coords_random, combined = TRUE)

log_cd_random <- matrix(temp%*%chol(S_random) + rep(c(1,2), each = n_points), ncol = 2)

informed_test_random <- fit_biwm(log_cd_random, coords_random, true_param, nus_vec)


generic_test_random <- fit_biwm(log_cd_random, coords_random, initial_guess_1, nus_vec)



exp_res_random <- 
  as.data.frame(
    rbind(
      true_param = round(true_param,5),
      initial_guess = round(initial_guess_1,5),
      estimated_generic = round(generic_test_random$theta,5),
      estimated_informed = round(informed_test_random$theta,5)
    )
  )


