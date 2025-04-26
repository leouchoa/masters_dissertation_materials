#' Guilherme este código/documento serve para ajudar a entender o problema principal do momento. A ideia é simular um simular um PG pelo `RandomFields` com RMbiwm em que $\sigma^2_1 = 1, \sigma^2_2 = 1.5, a = 1, \rho = 0.4, \nu_1 = \nu_2 = 1$.
#'
#' É importante citar também que $\nu_{red}$ é dado por $ν_{12} =ν_{21} = 0.5 (ν_{11} + ν_{22}) * ν_{red}$ e então para que tenhamos o modelo reduzido é preciso que $ν_{red} = 1$.
#'
#'
#'

# suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(BivMaternEstim))
suppressPackageStartupMessages(library(RandomFields))
# source("../../For_15-05/gradient/aux_fun/plot_rf_biwm.R")
set.seed(123)

# theme_set(theme_minimal())
# options(ggplot2.continuous.colour="viridis")

cov_struct <- RMbiwm(nudiag=c(1, 1), nured=1, rhored=0.4, cdiag=c(1, 1.5),s=c(1, 1, 1))

x_grid <- seq(-2, 2, 0.1)
load("../../For_29-05/gradient_problem/gp_and_idx.rda")
# gp <- RFsimulate(cov_struct, x_grid, x_grid)

# plot_rf_biwm(as.data.frame(gp),x_grid,x_grid)

coords <- as.data.frame(coordinates(gp))
obs_val <- as.data.frame(gp)

# idx <- sample(nrow(coords),100)

coords_subset <- coords[idx,]
obs_val_subset <- obs_val[idx,]

# ggplot(coords_subset,aes(coords.x1,coords.x2)) +
#   geom_point()

true_theta <- c(1,1.5,1,0.4)

theta_start_all_above <- c(1.5,3,4,0.7)

theta_start_all_bellow <- c(0.3,0.9,0.2,0.1)

theta_start_all_mixed_1 <- c(1.5,0.2,4,-0.3)


#' Comandos úteis pra ganhar tempo:
#'
#' - `debugonce(LLike_biwm);fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))`
#'
#' - `debugonce(BivMaternEstim:::block_LLike_biwm_grad);fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))`
#' 
#' - `block_LLike_biwm_grad(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))`
#'
#'
#' # Idéia do que eu tentei
#' 
#' 
#' 
#' Considerar os casos em que os valores iniciais estão todos acima do verdadeiro (`theta_start_all_above`), abaixo do verdadeiro (`theta_start_all_bellow`) ou misturados (`theta_start_all_mixed_1`). Com isso acompanhar $\theta$ e seu gradiente para tentar explicar o motivo do algorítmo estar com erro.
#' 
#' 
#' 
#' # Resumo
#' 
#' Não parece que a situação mudou muito depois da última reunião (06-06-2020), apesar do gráfico de contorno indicar que agora estamos otimizando a verossimilhança "correta".
#' 
#' ## Caso Todos acima
#' 
#' **Conclusão**:
#' 
#' 1. Me parece que o caso problemático aqui é o $\rho$ (fronteiras dele são +/- 1). Ele começou na direção errada e não mudou. Porém os outros também foram na direção da fronteira e não voltaram. Não dá pra afirmar que $\rho$ é o único problema.
#' 
#' Entretanto ele é bem problemático. Veja isto com o comando `res_test <- fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = c(theta_start_all_above[1:3],-0.2),nus = rep(1,2))`
#' 
#' 
#' 
#' 
#' **Resultado do otimizador:**
#' 
#' Comando usado: `res_all_above <- fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start_all_above,nus = rep(1,2))`
#' 
#' 
#' 
#' theta  1.5 3 4 0.7 
#' 
#' gradient  -20.68234 -12.20152 -36.8453 47.02203 
#' 
#' theta  1.21322 2.42626 3.23495 0.7572 
#' 
#' gradient  -22.58068 -13.92903 -51.46208 38.92481 
#' 
#' theta  0.06612 0.13129 0.17473 0.98601 
#' 
#' Error in optim(theta0, fn = LLike_biwm, gr = block_LLike_biwm_grad, method = "L-BFGS-B",  : 
#'                  L-BFGS-B needs finite values of 'fn'
#' 
#' 
#' ## Caso Todos Abaixo
#' 
#' **Conclusão:**
#' 
#' 1. Eles começam na direção e talvez isso seja fatal. Aqui o problema é o parâmetro $a$: por ser muito grande ele é o responsável por chegar à fronteira primeiro. Isto pode ser conferido com o comando `res_all_bellow <- fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = c(theta_start_all_below[1:3],-0.2),nus = rep(1,2))`
#' 
#' 2. Notar que a direção inicial foi exatamente a mesma do caso anterior
#' 
#' **Resultado do otimizador:**
#' 
#' Comando usado: `res_all_bellow <- fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start_all_bellow,nus = rep(1,2))`
#' 
#' 
#' theta  0.3 0.9 0.2 0.1 
#' 
#' gradient  -5.00702 -33.64344 -609.9887 3.39877 
#' 
#' theta  0.07368 0.21953 0.04937 0.78047 
#' 
#' Error in optim(theta0, fn = LLike_biwm, gr = block_LLike_biwm_grad, method = "L-BFGS-B",  : 
#'                  L-BFGS-B needs finite values of 'fn'
#' 
#' 
#' 
#' ## Caso Tudo Misturado
#' 
#' **Conclusão:**
#' 
#' 1. Curiosamente aqui ele acerta as direções e os primeiros passos são bons, mas depois ele volta pra onde estava e não anda pra lugar algum.
#' 
#' 2. Outro fato interessante é que o valor do gradiente não está muito longe daqueles encontrados nos outros testes. Mas aqui a atualização não ocorre, enquanto que nos outros casos eles vão pra fronteira muito rapidamente.
#' 
#' 3. 
#' 
#' 4. 
#' 
#' 
#' Comando usado: `res_all_mixed <- fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start_all_mixed_1,nus = rep(1,2))`
#' 
#' 
#' `r knitr::kable(as.table(rbind(true_theta,theta_start = theta_start_all_mixed_1,theta_hat = theta_start_all_mixed_1)))`
#' 
#' **Resultados das primeiras iterações:**
#' 
#' theta  1.5 0.2 4 -0.3 
#' 
#' gradient  -26.71931 290.6845 -54.36438 41.90233 
#' 
#' theta  1.49484 1.19988 3.98624 -0.29553 
#' 
#' gradient  -26.65159 -26.82395 -35.00898 -15.46164 
#' 
#' theta  1.49875 0.44333 3.99665 -0.29891 
#' 
#' gradient  -26.64757 -3.12186 -41.50066 3.70776 
#' 
#' theta  1.49974 0.25022 3.99931 -0.29978 
#' 
#' gradient  -26.69307 145.3303 -49.64324 27.86688 
#' 
#' theta  1.49996 0.20853 3.99988 -0.29996 
#' 
#' gradient  -26.71406 257.5295 -53.40168 39.03916 
#' 
#' theta  1.49999 0.2013 3.99998 -0.29999 
#' 
#' gradient  -26.71849 285.3268 -54.21229 41.44995
#' 
#' 
#' **Duas últimas iterações**
#' 
#' theta  1.5 0.2 4 -0.3 
#' 
#' gradient  -26.71931 290.6845 -54.36438 41.90233 
#' 
#' theta  1.5 0.2 4 -0.3 
#' 
#' gradient  -26.71931 290.6845 -54.36438 41.90233 
#' 
#' theta  1.5 0.2 4 -0.3 
#' 
#' gradient  -26.71931 290.6845 -54.36438 41.90233 
#' 
#' 
#' **O algorítmo "convergiu" para (1.5,0.2,4.0,-0.3).**
#' 
