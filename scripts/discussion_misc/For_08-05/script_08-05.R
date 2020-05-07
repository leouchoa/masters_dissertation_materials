library(ggplot2)
library(BivMaternEstim)
library(RandomFields)
library(gstat) #cross_variog

theme_set(theme_minimal())
options(ggplot2.continuous.colour="viridis")

# ----- Random Fields Simulation and Estimation ----
# takehome msg: the algorithm updates are too tiny to take a meaningfull step

set.seed(123)
# RFoptions(seed=1) 
cov_struct <- RMbiwm(nudiag=c(1, 1), nured=1, rhored=0.5, cdiag=c(1, 1.5),s=c(1, 1, 1)); plot(cov_struct)

x_grid <- seq(-2, 2, 0.1)
gp <- RFsimulate(cov_struct, x_grid, x_grid)

X11();plot(gp)

coords <- as.data.frame(coordinates(gp))
obs_val <- as.data.frame(gp)

idx <- sample(nrow(coords),100)

coords_subset <- coords[idx,]
obs_val_subset <- obs_val[idx,]

ggplot(coords_subset,aes(coords.x1,coords.x2)) + 
  geom_point()

theta_start <- c(1.5,3,4,0.5)

res <- fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))

debugonce(LLike_biwm);fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))

# O tamanho do passo tá minúsculo e ele não tá saindo de onde está.  Para descobrir isto basta dar `debugonce(LLike_biwm)` uma ou duas vezes dentro do algoritmo anterior. Dá prá notar que as atualizacoes sao minimas   

debugonce(BivMaternEstim:::block_LLike_biwm_grad);fit_biwm(obs_matrix = obs_val_subset,coords_matrix = obs_val_subset,theta0 = theta_start,nus = rep(1,2))

# ----- Original Dataset Estimationm ----
# takehome msgs:
# 
# - can't invert block cov matrices. I guess it's because the same coordinates are getting in the way
# 
# - imho the cross-variogram says we there's no meaning in using cross-covariance

dados <- na.omit(read.csv("../data/alr_transformed_with_locations.csv"))

coordinates(dados) <-  ~coord_x+coord_y
g <- gstat(NULL, 
           id = "alr_comp1", 
           formula = alr_comp1 ~ 1,
           data = dados)
g <- gstat(g, 
           id = "alr_comp2", 
           formula = alr_comp2 ~ 1,
           data = dados)

fct_labels <- c(alr_comp1.alr_comp2 = "Variograma Cruzado",
                alr_comp1 = "Variograma Direto alr 1",
                alr_comp2 = "Variograma Direto alr 2")

cross_vario <- variogram(g)
ggplot(cross_vario[,c(2,3,6)], aes(dist,gamma)) + 
  geom_point() + 
  facet_wrap(~id, 
             labeller = labeller(id = fct_labels))


theta_st <- c(var(dados$alr_comp1)
                 ,var(dados$alr_comp2)
                 , 0.3, 
                 cor(dados$alr_comp1,dados$alr_comp2))

res_2 <- fit_biwm(obs_matrix = dados[,1:2],coords_matrix = dados[,3:4],theta0 = theta_st,nus = c(0.5,1))


# ----- Original Dataset Estimationm ----
# takehome msgs:
# 
# - imho the cross-variogram says we there's no meaning in using cross-covariance
# 
# - Can estimate but the algorithm updates are too tiny to take a meaningfull step

dados_unique <- na.omit(read.csv("../data/alr_transformed_with_locations_UNIQUE.csv"))

coordinates(dados_unique) <-  ~coord_x+coord_y
vgm_alr1 <- variogram(alr_comp1 ~ 1, data = dados_unique,cloud=FALSE)
vgm_alr2 <- variogram(alr_comp2 ~ 1, data = dados_unique,cloud=FALSE)

# gridExtra::grid.arrange(
#   plot(vgm_alr1, main = "Variograma para log(areia/silte)"),
#   plot(vgm_alr2, main = "Variograma para log(argila/silte)"),
#   ncol = 2)

g_unique <- gstat(NULL, 
           id = "alr_comp1", 
           formula = alr_comp1 ~ 1,
           data = dados_unique)
g_unique <- gstat(g, 
           id = "alr_comp2", 
           formula = alr_comp2 ~ 1,
           data = dados_unique)

cross_vario_unique <- variogram(g)
ggplot(cross_vario_unique[,c(2,3,6)], aes(dist,gamma)) + 
  geom_point() + 
  facet_wrap(~id, 
             labeller = labeller(id = fct_labels))
# plot(cross_vario, pl=T)

# Sanity check

gridExtra::grid.arrange(
  ggplot(cross_vario[,c(2,3,6)], aes(dist,gamma)) + 
    geom_point() + 
    facet_wrap(~id, 
               labeller = labeller(id = fct_labels)),
  
  ggplot(cross_vario[,c(2,3,6)], aes(dist,gamma)) + 
    geom_point() + 
    facet_wrap(~id, 
               labeller = labeller(id = fct_labels))
  
)
