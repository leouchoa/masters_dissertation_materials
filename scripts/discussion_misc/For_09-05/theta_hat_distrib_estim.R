library(BivMaternEstim)
set.seed(1234)
n_points <- 7
nus_vec <- c(0.5, 0.5)

true_param <- c(1.5,0.7,1.3,0.4)
true_param <- c(1.5,0.7,1.3,-0.4)
true_param <- c(1.5,0.7,1.3,0)
initial_guess_1 <-  c(0.5,1.2,4,-0.2)
initial_guess_2 <-  c(0.5,1.2,4,0.4) #true rho


coords_random <- matrix(runif(2 * (n_points^2) ,-1,1), ncol = 2)

S_random <- sigma_assembler_biwm(sigmas = true_param[1:2], a = true_param[3], rho = true_param[4], nus = nus_vec, coords_matrix = coords_random, combined = TRUE)

S_chol <- chol(S_random)

get_distrib <- function(starting_point,n_replicates){
  
  
  res_generic_vec_01 <- replicate(n_replicates,{
    
    temp <- rnorm(2* (n_points^2))
    
    
    log_cd_random <- matrix(temp %*% S_chol + rep(c(1,2), each = n_points), ncol = 2)
    
    fit_biwm(log_cd_random, coords_random, starting_point, nus_vec)$theta
    
  })
  
  
  
  par(mfrow = c(2,2), pty = "s")
  
  hist(res_generic_vec_01[1,],xlab = "");abline(v = true_param[1], col = "red",lwd = 5)
  hist(res_generic_vec_01[2,],xlab = "");abline(v = true_param[2], col = "red",lwd = 5)
  hist(log(res_generic_vec_01[3,]),xlab = "");abline(v = log(true_param[3]), col = "red",lwd = 5)
  hist(-1 * res_generic_vec_01[4,],xlab = "");abline(v = true_param[4], col = "red",lwd = 5)
  
  return(res_generic_vec_01)
}


resultado <- get_distrib(initial_guess_2,30)


