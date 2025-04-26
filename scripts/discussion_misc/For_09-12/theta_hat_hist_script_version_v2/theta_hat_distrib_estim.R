library(BivMaternEstim)
set.seed(1234)
n_points <- 20
nus_vec <- c(0.5, 0.5)

temp <- rnorm(2* (n_points))

true_param <- c(1.5,0.7,1.3,-0.4)
initial_guess <-  c(0.5,5,4,0.2)

n_replicates <- 30

res_generic_vec_01 <- replicate(n_replicates,{
  
  coords_random <- matrix(runif(2 * (n_points) ,-1,1), ncol = 2)
  
  S_random <- sigma_assembler_biwm(sigmas = true_param[1:2], a = true_param[3], rho = true_param[4], nus = nus_vec, coords_matrix = coords_random, combined = TRUE)
  
  S_chol <- chol(S_random)
  
  log_cd_random <- matrix(temp %*% S_chol + rep(c(1,5), each = n_points), ncol = 2)
  
  fit_biwm(log_cd_random, coords_random, initial_guess, nus_vec)$theta
  
})



par(mfrow = c(2,2), pty = "s")

hist(res_generic_vec_01[1,],xlab = "");abline(v = true_param[1], col = "red",lwd = 5)
hist(res_generic_vec_01[2,],xlab = "");abline(v = true_param[2], col = "red",lwd = 5)
hist(res_generic_vec_01[3,],xlab = "");abline(v = true_param[3], col = "red",lwd = 5)
hist(res_generic_vec_01[4,],xlab = "");abline(v = true_param[4], col = "red",lwd = 5)




