#' Guilherme este serve para ajudar a entender o problema principal do momento. A ideia é simular um simular um PG pelo RandomFields com RMbiwm em que $\sigma^2_1 = 1, \sigma^2_2 = 1.5, a = 1, \rho = 0.7, \nu_1 = \nu_2 = 1$. 
#' 
#' É importante citar também que $\nu_{red}$ é dado por $ν_{12} =ν_{21} = 0.5 (ν_{11} + ν_{22}) * ν_{red}$ e então para que tenhamos o modelo reduzido é preciso que $ν_{red} = 1$.
#' 
#' 
#' 

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(BivMaternEstim))
suppressPackageStartupMessages(library(RandomFields))
source("../aux_fun/plot_rf_biwm.R")

theme_set(theme_minimal())
options(ggplot2.continuous.colour="viridis")

set.seed(123)
cov_struct <- RMbiwm(nudiag=c(1, 1), nured=1, rhored=0.7, cdiag=c(1, 1.5),s=c(1, 1, 1))

x_grid <- seq(-2, 2, 0.1)
gp <- RFsimulate(cov_struct, x_grid, x_grid)

plot_rf_biwm(as.data.frame(gp),x_grid,x_grid)

coords <- as.data.frame(coordinates(gp))
obs_val <- as.data.frame(gp)

idx <- sample(nrow(coords),100)

coords_subset <- coords[idx,]
obs_val_subset <- obs_val[idx,]

ggplot(coords_subset,aes(coords.x1,coords.x2)) + 
  geom_point()

theta_start <- c(1.5,3,4,0.5)

res <- fit_biwm(obs_matrix = obs_val_subset,
                coords_matrix = obs_val_subset,
                theta0 = theta_start,
                nus = rep(1,2),
                verbosity = 6)

# debugonce(LLike_biwm);fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))
# 
#' O tamanho do passo tá minúsculo e ele não tá saindo de onde está. Para descobrir isto usei `debugonce(LLike_biwm)` duas vezes dentro do algoritmo anterior. Também é importante citar que com `debugonce(BivMaternEstim:::block_LLike_biwm_grad)` dá pra notar que os gradientes estão com valores muito altos.
# 
# debugonce(BivMaternEstim:::block_LLike_biwm_grad);fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))