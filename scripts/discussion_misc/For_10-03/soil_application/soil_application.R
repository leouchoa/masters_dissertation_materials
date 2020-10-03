library(BivMaternEstim)
library(gstat)
library(sp)
library(ggplot2)
library(tidyr)
library(ggmap)
library(patchwork)

alr_transformed_with_locations_UNIQUE <- read.csv("~/Documents/git/Mvt_Spat_Stat/scripts/discussion_misc/data/alr_transformed_with_locations_UNIQUE.csv", row.names=NULL)

soil_dts <- alr_transformed_with_locations_UNIQUE[,1:4]

# nrow(unique(soil_dts)) == nrow(soil_dts)

coordinates(soil_dts) <-  ~coord_x + coord_y

g <- gstat(NULL, 
           id = "alr_comp1", 
           formula = alr_comp1 ~ 1,
           data = soil_dts)
g <- gstat(g, 
           id = "alr_comp2", 
           formula = alr_comp2 ~ 1,
           data = soil_dts)

fct_labels <- c(alr_comp1.alr_comp2 = "Variograma Cruzado",
                alr_comp1 = "Variograma Processo 1",
                alr_comp2 = "Variograma Processo 2")

cross_vario <- variogram(g)
ggplot(cross_vario[,c(2,3,6)], aes(dist,gamma)) + 
  geom_point() + 
  facet_wrap(~id, 
             labeller = labeller(id = fct_labels)
  ) +
  theme_minimal()

initial_point <- c(1,1,0,0.5)
nus_vec <- c(1,1)

soil_dts_fit <- fit_biwm(
  obs_matrix = soil_dts@data,
  coords_matrix = coordinates(soil_dts),
  theta0 = initial_point,
  nus = nus_vec
  )

coords_range <- apply(coordinates(soil_dts),2,range)

krig_grid <- 
  expand.grid(
    seq(coords_range[1,1],coords_range[2,1],length.out = 30),
    seq(coords_range[1,2],coords_range[2,2],length.out = 30)
  )

soil_dts_preds <- compositional_biwm_krig(
  biwm_fit = soil_dts_fit,
  krig_locations = krig_grid,
  fit_locations = coordinates(soil_dts),
  obs_matrix = soil_dts@data,
  nus = nus_vec
  )

soil_dts_preds_plot <- soil_dts_preds %>% 
  gather("key","Valor",1:3) %>% 
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller = 
               labeller(key =c(comp_01 = "Composição 1",
                               comp_02 = "Composição 2",
                               comp_03 = "Composição 3")),
             ncol = 1) + 
  labs(x = "", y = "") +
  scale_fill_viridis_c()

soil_dts_preds_plot

# i think it is as shitty as it is bcs of the nugget

ctb_map_image <- readRDS("ctb_map_image_file")

soil_dts_map_plot <- ggmap(ctb_map_image) + 
  geom_point(data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")), 
             mapping = aes(x = coord_x,y = coord_y))

soil_dts_preds_plot + soil_dts_map_plot

# --------- New data points ----------


NEW_alr_transformed_with_locations_UNIQUE <- read.csv("~/Documents/git/Mvt_Spat_Stat/scripts/discussion_misc/data/NEW_alr_transformed_with_locations_UNIQUE.csv", row.names=NULL)[,c(4:8,10,11)]

NEW_obs_soil_dts_preds <- compositional_biwm_krig(
  biwm_fit = soil_dts_fit,
  krig_locations = NEW_alr_transformed_with_locations_UNIQUE[,6:7],
  fit_locations = coordinates(soil_dts),
  obs_matrix = soil_dts@data,
  nus = nus_vec
)

# the error is large
colSums(NEW_obs_soil_dts_preds[,1:3] - NEW_alr_transformed_with_locations_UNIQUE[,1:3])

apply(
  X = NEW_obs_soil_dts_preds[,1:3] - NEW_alr_transformed_with_locations_UNIQUE[,1:3],
  MARGIN = 2,
  FUN = hist
  )