#' Gráficos de Contorno para a Log-verossimilhança, dados valores de parâmetros, que são $\sigma^2_1,\sigma^2_2,a,\rho$. Então para fazer o gráfico basta construir uma malha com a função `construct_grid` e chamar a função `plot_llike_contour`, cujos dois primeiros argumentos devem ter um dos nomes `sigma2_1,sigma2_1,a,rho` e o terceiro argumento deve ser a malha. 
#' 
#' Exemplos :
#' 
#' 

library(BivMaternEstim)
library(ggplot2)
source("contour_funs.R")
set.seed(3)
options(ggplot2.continuous.colour="viridis")


n <- 40
coords <- matrix(runif(2*n), ncol = 2)
temp <- rnorm(2*n)
nus_vec <- c(0.5,0.5)


grid_01 <- construct_grid(
  seq(0.01,10,length.out = 15),
  1,
  seq(0.01,10,length.out = 15),
  0.2)

plot_llike_contour("sigma2_1","a",grid_01)

grid_02 <- construct_grid(
  1,
  1,
  seq(0.01,10,length.out = 15),
  seq(0,0.7,length.out = 15)
)

plot_llike_contour("a","rho",grid_02)

grid_03 <- construct_grid(
  seq(0.01,10,length.out = 15),
  seq(0.01,10,length.out = 15),
  1,
  0.2)

plot_llike_contour("sigma2_1","sigma2_2",grid_03)

grid_04 <- construct_grid(
  seq(0.01,10,length.out = 15),
  1,
  1,
  seq(-0.8,0.8,length.out = 15))

plot_llike_contour("sigma2_1","rho",grid_04)

grid_04 <- construct_grid(
  seq(0.01,10,length.out = 15),
  1,
  1,
  seq(-0.8,0.8,length.out = 15))

plot_llike_contour("sigma2_1","rho",grid_04)

