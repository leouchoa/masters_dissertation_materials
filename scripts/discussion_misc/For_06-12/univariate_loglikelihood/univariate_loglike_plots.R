#' Gráficos Univariados da Log-Verossimilhança
#' 
#' Guilherme aqui fiz isso com o intuito de ver melhor como está o formato da verossimilhança e o quão deslocado ela está em cada parâmetro, pois achei que nos gráficos de contorno não muito claro. 
#' 
#' Pra isso eu simulei um Processo Gaussiano com $\sigma^2_1 = 1, \sigma^2_2 = 1, a = 2, \rho = 0.5, \nu_1 = \nu_2 = 0.5$ da mesma forma que é simulado nos gráficos de contorno.
#' 
#' 
#' 
#' ## Simulação do Processo

suppressPackageStartupMessages(library(BivMaternEstim))
source("univariate_loglike_funs.R")

n <- 100
coords <- matrix(runif(2*n), ncol = 2)
temp <- rnorm(2*n)
nus_vec <- c(0.5,0.5)
S <- sigma_assembler_biwm(sigmas = c(1, 1), a = 2, rho = 0.5,
                          nus = c(0.5, 0.5), coords_matrix = coords, combined = TRUE)

log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)



#' ## Gráficos para $\sigma^2_1$
#' 
#' 

grid_sigma2_1_notzoomed <- construct_grid(
  seq(0.01,2.5,length.out = 100),
  1,
  1,
  0.2)

grid_sigma2_1_zoomed <- construct_grid(
  seq(0.75,1.25,length.out = 100),
  1,
  1,
  0.2)


p_sigma2_1_notzoomed <- plot_univariate_llike("sigma2_1",grid_sigma2_1_notzoomed,log_cd,1)

p_sigma2_1_zoomed <- plot_univariate_llike("sigma2_1",grid_sigma2_1_zoomed,log_cd,1)

plot_both_zoom_llike(p_sigma2_1_notzoomed,p_sigma2_1_zoomed,"sigma2_1",2)


#' ## Gráficos para $\sigma^2_2$
#' 
#' 
#' Parece idiota fazer isso, mas é só pra garantir
#' 
#' 



grid_sigma2_2_notzoomed <- construct_grid(
  1,
  seq(0.01,2.5,length.out = 100),
  1,
  0.2)

grid_sigma2_2_zoomed <- construct_grid(
  1,
  seq(0.75,1.25,length.out = 100),
  1,
  0.2)


p_sigma2_2_notzoomed <- plot_univariate_llike("sigma2_2",grid_sigma2_2_notzoomed,log_cd,1)

p_sigma2_2_zoomed <- plot_univariate_llike("sigma2_2",grid_sigma2_2_zoomed,log_cd,1)

plot_both_zoom_llike(p_sigma2_2_notzoomed,p_sigma2_2_zoomed,"sigma2_2",2)



#' ## Gráficos para $a$
#' 
#' 
#' 


grid_a_notzoomed <- construct_grid(
  1,
  1,
  seq(0.01,2.5,length.out = 100),
  0.2)

grid_a_zoomed <- construct_grid(
  1,
  1,
  seq(0.75,1.25,length.out = 100),
  0.2)


p_a_notzoomed <- plot_univariate_llike("a",grid_a_notzoomed,log_cd,2)

p_a_zoomed <- plot_univariate_llike("a",grid_a_zoomed,log_cd,2)

plot_both_zoom_llike(p_a_notzoomed,p_a_zoomed,"a",2)







#' ## Gráficos para $\rho$
#' 
#' 

grid_rho_notzoomed <- construct_grid(
  1,
  1,
  2,
  seq(-0.95,0.95,length.out = 200))

grid_rho_zoomed <- construct_grid(
  1,
  1,
  2,
  seq(-0.25,0.75,length.out = 100))

p_rho_notzoomed <- plot_univariate_llike("rho",grid_rho_notzoomed,log_cd,0.5)

p_rho_zoomed <- plot_univariate_llike("rho",grid_rho_zoomed,log_cd,0.5)

plot_both_zoom_llike(p_rho_notzoomed,p_rho_zoomed,"rho",2)