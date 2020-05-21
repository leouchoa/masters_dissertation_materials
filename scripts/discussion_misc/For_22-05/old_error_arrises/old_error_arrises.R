#' Guilherme este serve para ajudar a entender o problema principal do momento. A ideia é simular um simular um PG pelo RandomFields com RMbiwm em que $\sigma^2_1 = 1, \sigma^2_2 = 1.5, a = 1, \rho = 0.7, \nu_1 = \nu_2 = 1$. 
#' 
#' É importante citar também que $\nu_{red}$ é dado por $ν_{12} =ν_{21} = 0.5 (ν_{11} + ν_{22}) * ν_{red}$ e então para que tenhamos o modelo reduzido é preciso que $ν_{red} = 1$.
#' 
#' 
#' 

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(BivMaternEstim))
suppressPackageStartupMessages(library(RandomFields))
source("../../For_15-05/gradient/aux_fun/plot_rf_biwm.R")
set.seed(123)

theme_set(theme_minimal())
options(ggplot2.continuous.colour="viridis")

cov_struct <- RMbiwm(nudiag=c(1, 1), nured=1, rhored=0.7, cdiag=c(1, 1.5),s=c(1, 1, 1))

x_grid <- seq(-2, 2, 0.1)
gp <- RFsimulate(cov_struct, x_grid, x_grid)

plot_rf_biwm(as.data.frame(gp),x_grid,x_grid)

coords <- as.data.frame(coordinates(gp))
obs_val <- as.data.frame(gp)

idx <- sample(nrow(coords),100)

coords_subset <- coords[idx,]
obs_val_subset <- obs_val[idx,]

# ggplot(coords_subset,aes(coords.x1,coords.x2)) + 
#   geom_point()

theta_start <- c(1.5,3,4,0.5)


# Comentei o código abaixo só para compilar o documento.
#
# res <- fit_biwm(obs_matrix = obs_val_subset,
#                 coords_matrix = obs_val_subset,
#                 theta0 = theta_start,
#                 nus = rep(1,2),
#                 verbosity = 6)

#' Comandos úteis pra ganhar tempo:
#' 
#' - `debugonce(LLike_biwm);fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))`
#'
#' - `debugonce(BivMaternEstim:::block_LLike_biwm_grad);fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))`
#'
#' # Resumo
#' 
#' Consegui corrigir um problema no cálculo do gradiente e a situação mudou. Agora o problema volta a ser o de anteriormente em que a função objetivo tá explodindo na hora de calcular $\log(|\boldsymbol\Sigma_{\boldsymbol \theta}|)$. O erro é
#' 
#' `Error in optim(theta0, fn = LLike_biwm, gr = block_LLike_biwm_grad, method = "L-BFGS-B",  :  L-BFGS-B needs finite values of 'fn'`
#' 
#' # O que eu já tentei
#' 
#' 1. Conferir as inversas calculada por blocos e com o `solve` e o resultado é que elas são quase iguais. Dentro da `debugonce(LLike_biwm)` eu fiz
#' 
#' - `autocov_test <- sigma_assembler_biwm(sigmas = sigmas,a = as,rho = rho,nus = nus,coords_matrix = coords_matrix,combined = TRUE)`
#' - `autocov_test_inv <- solve(autocov_test)`
#' - `image(autocov_test_inv - block_inv(autocov_matrix,combined = TRUE))`
#' 
#' 
#' 
#' 2. Comparar o $\log(|\boldsymbol\Sigma_{\boldsymbol \theta}|)$ por blocos usando `block_log_det()` e com `log(det())`. O resultado foi que elas dão o mesmo resultado. Pra chegar nessa conclusão eu aproveitei as matrizes que foram criadas no tópico 1. dentro de `debugonce(LLike_biwm)` e fiz
#' 
#' - `round(block_log_det(autocov_matrix) - log(det(autocov_test)), 10)`
#' 
#' Na terceira vez que usei `debugonce(Llike_biwm)` o log do determinante explodiu pra `-Inf` tanto por `block_log_det(autocov_matrix)` como por `log(det(autocov_test))`.
#' 
#' 
#' 
#' 3. Acompanhar o gradiente da função nas iterações com `debugonce(BivMaternEstim:::block_LLike_biwm_grad)`. Na segunda iteração o gradiente não havia explodido. Nestas duas corridas ambos os gradientes de $\sigma^2_1$ e $\sigma^2_2$ concordavam em sinal e na ordem de magnitude.

