library(BivMaternEstim)
library(ggplot2)
library(gridExtra)
theme_set(theme_minimal())


set.seed(1234)

true_param <- c(1.5,0.7,1.3,-0.4)
initial_guess <-  c(0.5,5,4,0.2)
n_points <- 80
nus_vec <- c(0.5,0.5)
n_replicates <- 30


coords_random <- matrix(runif(2 * (n_points) ,-1,1), ncol = 2)

S_random <- sigma_assembler_biwm(sigmas = true_param[1:2], 
                                 a = true_param[3], 
                                 rho = true_param[4], 
                                 nus = nus_vec, 
                                 coords_matrix = coords_random, 
                                 combined = TRUE
)

S_chol <- chol(S_random)

res_generic_vec_01 <- replicate(n_replicates,{
  
  temp <- rnorm(2* (n_points))
  
  
  log_cd_random <- matrix(temp %*% S_chol + rep(c(1,2), each = n_points), ncol = 2)
  
  fit_biwm(log_cd_random, coords_random, initial_guess, nus_vec)$theta
  
})


res_generic_vec_01 <- 
  setNames(
    as.data.frame(
      t(
        res_generic_vec_01
      )
    ),
    c("sigma2_1","sigma2_2","a","rho")
  )

p_sigma2_1 <- ggplot(res_generic_vec_01,aes(sigma2_1)) + 
  geom_histogram() + 
  labs(
    x = expression(sigma[1]^2),
    y = "Contagem",
    title = expression("Histograma de" ~ sigma[1]^2 ~ "Estimado")
  )  + 
  geom_vline(xintercept = true_param[1])


p_sigma2_2 <- ggplot(res_generic_vec_01,aes(sigma2_2)) + 
  geom_histogram() + 
  labs(
    x = expression(sigma[2]^2),
    y = "Contagem",
    title = expression("Histograma de" ~ sigma[2]^2 ~ "Estimado")
  )   + 
  geom_vline(xintercept = true_param[2])

p_rho <- ggplot(res_generic_vec_01,aes(a)) + 
  geom_histogram() + 
  labs(
    x = expression(a),
    y = "Contagem",
    title = expression("Histograma de" ~ a~ "Estimado")
  )   + 
  geom_vline(xintercept = true_param[3])


p_a <- ggplot(res_generic_vec_01,aes(rho)) + 
  geom_histogram() + 
  labs(
    x = expression(rho),
    y = "Contagem",
    title = expression("Histograma de" ~ rho~ "Estimado")
  )   + 
  geom_vline(xintercept = true_param[4])

grid.arrange(
  p_sigma2_1,
  p_sigma2_2,
  p_rho,
  p_a
)



