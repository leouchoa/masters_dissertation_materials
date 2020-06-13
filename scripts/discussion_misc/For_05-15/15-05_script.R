# ----- Gaussian Process ---

#' Guilherme corrigi o código anterior e comentei um pouco esse aqui pra eu não desperdiçar o seu tempo.
#' 
#' A idéia é eu gostaria que você desse uma olhada aqui nos resultados só pra conferir que eu não fiz mais cagada. Já basta semana passada e retrasada.
#' 
#' A imagem derivada de `source("biwm_plot.R")` mostra o exemplo da figura 04 que eu coloquei no overleaf. Dela eu refiz os exemplos e agora eles de fato seguem a estrutura que está lá. Esta estrutura de fato permite a formulação reduzida do modelo já que $|\rho| \leq \frac{\sqrt{(\nu_1 \nu_2)}}{(\nu_1 + \nu_2) / 2}$ pois 0.5 + 1.5 / 2 = 1 e $\sqrt{(0.5 *1.5)}$ = `r sqrt(0.5*1.5)`.
#'
#'Os gráficos posteriores são os resultados obtidos a partir da simulação para o processo baseado na tal figura 04. 
#'



source("biwm_plot.R")
library(RandomFields)
library(gstat)
library(sp)
library(tidyverse)
theme_set(theme_minimal())

set.seed(123)
dists <- seq(0,5,0.1)
biwm_plot(dists,1,3,1,0.8,c(0.5,1.5,1.0))# diss biwm_covs-struct example

# x <- y <- seq(-10, 10, 0.5)
# model <- RMbiwm(nudiag=c(0.5, 1.5), nured=1, rhored=0.8, cdiag=c(1, 3), s=c(1, 1, 1))
# 
# biwm_sim <- RFsimulate(model, x, y)
# 
# my_grid <- expand.grid(x,y)
# biwm_sim_df <- data.frame(coord_x = my_grid$Var1,coord_y = my_grid$Var2, process_1 = biwm_sim$variable1, process_2 = biwm_sim$variable2)

load("biwm_sim_&_idx.rda")


biwm_sim_df %>% 
  gather("key","Valor",3:4) %>% 
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller = 
               labeller(key =c(process_1 = "Processo 1",
                               process_2 = "Processo 2"))) + 
  labs(x = "", y = "") +
  scale_fill_viridis_c()

biwm_sim_df_comp <- cbind(biwm_sim_df[,1:2],compositions::alrInv(biwm_sim_df[,3:4])) %>% 
  setNames(c("coord_x",
             "coord_y",
             paste0("comp_",1:3)))

# ----- Gaussian Process Subset Estimation and Variogram ------

biwm_sim_df_comp %>% 
  gather("key","Valor",3:5) %>% 
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller = 
               labeller(key =c(comp_1 = "Composição 1",
                               comp_2 = "Composição 2",
                               comp_3 = "Composição 3"))) + 
  labs(x = "", y = "") +
  scale_fill_viridis_c()

# idx <- sample(nrow(biwm_sim_df),150)
# 
# save(biwm_sim_df,idx,file = "biwm_sim_&_idx.rda")

biwm_sim_df_sampled <- biwm_sim_df[idx,]

biwm_sim_df_sampled %>% 
  gather("key","val",-(coord_x:coord_y)) %>% 
  ggplot(aes(coord_x,coord_y,color = val)) + 
  geom_point() +
  facet_wrap(~key)

# --- Variogram ---

coordinates(biwm_sim_df_sampled) <-  ~coord_x+coord_y


g <- gstat(NULL, 
           id = "process_2", 
           formula = process_2 ~ 1,
           data = biwm_sim_df_sampled)

g <- gstat(g, 
           id = "process_1", 
           formula = process_1 ~ 1,
           data = biwm_sim_df_sampled)

fct_labels <- c(process_2.process_1 = "Variograma Cruzado",
                process_1 = "Variograma Direto alr 1",
                process_2 = "Variograma Direto alr 2")

cross_vario <- variogram(g)
ggplot(cross_vario[,c(2,3,6)], aes(dist,gamma)) + 
  geom_point() + 
  facet_wrap(~id, 
             labeller = labeller(id = fct_labels))

