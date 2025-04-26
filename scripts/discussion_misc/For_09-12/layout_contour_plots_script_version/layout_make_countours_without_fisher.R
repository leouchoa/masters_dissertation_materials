#' It's necessary to the script `contourField_funs.R` to be in the same folder as the file contains complementary functions
#' 
#' 
#' 
#' IT DOESNT WORK IF THE GRADIENT IS ALREADY PRE-MULTIPLIED BY `gamma_plot_val` because it will be really small. A **not** possible workaround its to pass `gamma_plot_val` with inverse scale 


suppressPackageStartupMessages(library(BivMaternEstim))
set.seed(123)
source("contourField_funs.R")
  


n =70
gamma_plot_val = 1e3
coords <- matrix(runif(2*n), ncol = 2)
temp <- rnorm(2*n)
nus_vec <- c(0.5,0.5)

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

layout(matrix(c(1, 2,
                3, 4), nr=2, byrow=T))
# hist(rnorm(25), col="VioletRed",main = "p1")
# hist(rnorm(25), col="VioletRed",main = "p2")
# hist(rnorm(25), col="VioletRed",main = "p3") 
# hist(rnorm(25), col="VioletRed",main = "p4")


grid_01 <- construct_grid(
  sigma2_1_grid,
  1,
  a_grid,
  0.2)



plot_llike_contour_field("sigma2_1","a",grid_df = grid_01,obs = log_cd,coords = coords,gamma_plot = gamma_plot_val)


grid_02 <- construct_grid(
  1,
  1,
  a_grid,
  rho_grid
)

plot_llike_contour_field("a","rho",grid_df = grid_02,obs = log_cd,coords = coords, gamma_plot = gamma_plot_val)






grid_03 <- construct_grid(
  sigma2_1_grid,
  sigma2_2_grid,
  2,
  0.5)


plot_llike_contour_field("sigma2_1","sigma2_2",grid_df = grid_03,obs = log_cd,coords = coords,gamma_plot = gamma_plot_val)




grid_04 <- construct_grid(
  sigma2_1_grid,
  1,
  1,
  rho_grid
)

plot_llike_contour_field("sigma2_1","rho",grid_df = grid_04,obs = log_cd,coords = coords,gamma_plot = gamma_plot_val)



