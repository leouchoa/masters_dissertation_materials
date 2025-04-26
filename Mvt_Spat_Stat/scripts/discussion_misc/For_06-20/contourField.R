#' Gráficos de Contorno para a Log-verossimilhança, dados valores de parâmetros, que são $\sigma^2_1,\sigma^2_2,a,\rho$. Então para fazer o gráfico basta construir uma malha com a função `construct_grid` e chamar a função `plot_llike_contour`, cujos dois primeiros argumentos devem ter um dos nomes `sigma2_1,sigma2_1,a,rho` e o terceiro argumento deve ser a malha. 
#' 
#' 
#' 

suppressPackageStartupMessages(library(BivMaternEstim))
# library(ggplot2)
source("contourField_funs.R")
set.seed(123)
# options(ggplot2.continuous.colour="viridis")
# theme_set(theme_minimal())


n <- 100
coords <- matrix(runif(2*n), ncol = 2)
temp <- rnorm(2*n)
nus_vec <- c(0.5,0.5)
S <- sigma_assembler_biwm(sigmas = c(1, 1), a = 2, rho = 0.5,
                          nus = c(0.5, 0.5), coords_matrix = coords, combined = TRUE)

log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)


grid_01 <- construct_grid(
  seq(0.8,2.5,length.out = 15),
  1,
  seq(0.8,2.5,length.out = 15),
  0.2)


plot_llike_contour(var_1 = "sigma2_1",
                   var_2 = "a",
                   grid_df = grid_01,
                   obs = log_cd,
                   true_param = c(1,1))

grid_02 <- construct_grid(
  1,
  1,
  seq(0.5,3,length.out = 15),
  seq(0.4,0.9,length.out = 15)
)

plot_llike_contour(var_1 = "a",
                   var_2 = "rho",
                   grid_df = grid_02,
                   obs = log_cd,
                   true_param = c(1,0.5))

grid_03 <- construct_grid(
  seq(0.25,1.5,length.out = 15),
  seq(0.25,1.5,length.out = 15),
  2,
  0.5)

plot_llike_contour(var_1 = "sigma2_1",
                   var_2 = "sigma2_2",
                   grid_df = grid_03,
                   obs = log_cd,
                   true_param = c(1,1))

grid_04 <- construct_grid(
  seq(0.25,1.5,length.out = 15),
  1,
  1,
  seq(0.4,0.8,length.out = 15))


plot_llike_contour(var_1 = "sigma2_1",var_2 = "rho",grid_df = grid_04,obs = log_cd,true_param = c(1,0.5))

