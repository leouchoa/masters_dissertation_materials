library(BivMaternEstim)

library(ggplot2)
library(patchwork)

n <- 70
coords <- matrix(runif(2*n), ncol = 2)
temp <- rnorm(2*n)
nus_vec <- c(0.5,0.5)
sample.int(123123,1)

true_theta = c(1,1,2,0.5)
sigma2_1_grid <- seq(0.5,2,length.out = 15)
sigma2_2_grid <- seq(0.5,2,length.out = 15)
a_grid <- seq(0.5,3,length.out = 15)
rho_grid <- seq(0.25,0.75,length.out = 15)

S <- sigma_assembler_biwm(sigmas = true_theta[1:2],
                          a = true_theta[3],
                          rho = true_theta[4],
                          nus = c(0.5, 0.5),
                          coords_matrix = coords,
                          combined = TRUE)

log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)

grid_01 <- construct_grid(
  sigma2_1_grid,
  1,
  a_grid,
  0.2)

grid_02 <- construct_grid(
  1,
  1,
  a_grid,
  rho_grid
)


grid_03 <- construct_grid(
  sigma2_1_grid,
  sigma2_2_grid,
  2,
  0.5)

grid_04 <- construct_grid(
  sigma2_1_grid,
  1,
  1,
  rho_grid
)



contour_sigma2_1_a <-
  plot_loglike_contour_ggplot("sigma2_1","a",grid_df = grid_01,obs = log_cd,coords = coords,nbins = 30,gamma_plot =1e3,true_param = true_theta[c(1,3)]) + theme(legend.position="none") +
  labs(
    x = expression(sigma[1]^2)
  )


contour_a_rho <- plot_loglike_contour_ggplot("a","rho",grid_df = grid_02,obs = log_cd,coords = coords,nbins = 30,gamma_plot =1e3,true_param = true_theta[c(3,4)]) +
  theme(legend.position="none") +
  labs(
    y = expression(rho)
  )



contour_sigma2_1_sigma2_2 <- plot_loglike_contour_ggplot("sigma2_1","sigma2_2",grid_df = grid_03,obs = log_cd,coords = coords,nbins = 30,gamma_plot =1e3,true_param = true_theta[c(1,2)]) +
  theme(legend.position="none") +
  labs(
    x = expression(sigma[1]^2),
    y = expression(sigma[2]^2)
  )



contour_sigma2_1_rho <- plot_loglike_contour_ggplot("sigma2_1","rho",grid_df = grid_04,obs = log_cd,coords = coords,nbins = 30,gamma_plot =1e3,true_param = true_theta[c(1,4)])  +
  theme(legend.position="none") +
  labs(
    x = expression(sigma[1]^2),
    y = expression(rho)
  )



contour_sigma2_1_a +
  contour_a_rho +
  contour_sigma2_1_sigma2_2 +
  contour_sigma2_1_rho
