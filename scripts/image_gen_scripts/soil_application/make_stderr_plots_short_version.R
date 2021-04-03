library(BivMaternEstim)
library(gstat)
library(sp)
library(ggplot2)
library(tidyr)
library(ggmap)
library(viridis)
library(soiltexture)
library(gridExtra)

#krig_grid of leng 150 is good

KRIG_GRID_LEN <- 95

alr_transformed_with_locations_UNIQUE <- read.csv("data/alr_transformed_with_locations_UNIQUE.csv", row.names=NULL)

soil_dts <- alr_transformed_with_locations_UNIQUE[,1:4]

# nrow(unique(soil_dts)) == nrow(soil_dts)

coordinates(soil_dts) <-  ~coord_x + coord_y
initial_point <- c(1,1,0.5,0)
nus_vec <- c(1,1)
nug_vec <- c(0.5,0.5)

cat("\nEstimating parameters........")

soil_dts_fit <- fit_biwm(
  obs_matrix = soil_dts@data,
  coords_matrix = coordinates(soil_dts),
  theta0 = initial_point,
  nus = nus_vec,
  nug_vec = nug_vec
)

cat("\nEstimation completed!")

coords_range <- apply(coordinates(soil_dts),2,range)

krig_grid <-
  expand.grid(
    seq(coords_range[1,1] - 0.05,coords_range[2,1] + 0.09,length.out = KRIG_GRID_LEN),
    seq(coords_range[1,2] - 0.05,coords_range[2,2] + 0.09,length.out = KRIG_GRID_LEN)
  )

cat("\nCreating compositional prediction map")

soil_dts_preds <- compositional_biwm_krig(
  biwm_fit = soil_dts_fit,
  krig_locations = krig_grid,
  fit_locations = coordinates(soil_dts),
  obs_matrix = soil_dts@data,
  nus = nus_vec,
  nug_vec
)

cat("\nCompositional prediction map created!")

names(soil_dts_preds) <- c("areia","argila","silte","coord_x","coord_y")

ctb_map_image <- readRDS("ctb_map_image_file")
ctb_map_image_unique_loc_constrained <- readRDS("ctb_map_image_unique_loc_constrained.rda")

ctb_map_image_unique_loc_constrained_v2 <- readRDS("ctb_map_image_unique_loc_constrained_v2.rda")


soil_dts_map_plot_v2 <- ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")), 
             mapping = aes(x = coord_x,y = coord_y))

krig_vals <- compositional_biwm_krig(
  biwm_fit = soil_dts_fit,
  krig_locations = krig_grid,
  fit_locations = coordinates(soil_dts),
  obs_matrix = soil_dts@data,
  nus = nus_vec,
  nug_vec,
  return_comp = FALSE
)

cat("\nComputing Approximate Covariance Inverse.....")

sigbar_hat_approx <- 
  get_approx_alrInv_cond_covvar(krig_values = krig_vals,loc_new = krig_grid,loc_obs = soil_dts@coords,biwm_fit = soil_dts_fit,nus = nus_vec,nug_vec = nug_vec)

cat("\nApproximate Covariance Inverse computed!!")


areia_std_err_dm <- ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(
    data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")),
    mapping = aes(x = coord_x,
                  y = coord_y,
                  size = alr_transformed_with_locations_UNIQUE$prop1_areia
    )
  ) +
  geom_tile(data = sigbar_hat_approx[,c(1,4,5)], 
            aes(coord_x,coord_y, fill = comp_01_sd, alpha = 0.2)
  ) + 
  scale_fill_viridis_c() +
  # scale_fill_gradientn(colours=viridis(7),na.value = "transparent",
  #                      breaks=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5),
  #                      limits=c(0,1))+
  # guides(alpha = FALSE, 
  #        size = FALSE) + 
  guides(alpha = FALSE) + 
  labs(
    x = "",
    y = "",
    title = "Desvio Padrão da Predição - Areia",
    size = "Areia Observada",
    fill = "Erro Padrão"
  ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

argila_std_err_dm <- ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(
    data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")),
    mapping = aes(x = coord_x,
                  y = coord_y,
                  size = alr_transformed_with_locations_UNIQUE$prop2_argila
    )
  ) +
  geom_tile(data = sigbar_hat_approx[,c(2,4,5)], 
            aes(coord_x,coord_y, fill = comp_02_sd, alpha = 0.2)
  ) + 
  scale_fill_viridis_c() +
  # scale_fill_gradientn(colours=viridis(7),na.value = "transparent",
  #                      breaks=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5),
  #                      limits=c(0,1))+
  # guides(alpha = FALSE, 
  #        size = FALSE) + 
  guides(alpha = FALSE) + 
  labs(
    x = "",
    y = "",
    title = "Desvio Padrão da Predição - Argila",
    size = "Argila Observada",
    fill = "Erro Padrão"
  ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

silte_std_err_dm <- ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(
    data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")),
    mapping = aes(x = coord_x,
                  y = coord_y,
                  size = alr_transformed_with_locations_UNIQUE$prop3_silte
    )
  ) +
  geom_tile(data = sigbar_hat_approx[,c(3,4,5)], 
            aes(coord_x,coord_y, fill = comp_03_sd, alpha = 0.2)
  ) + 
  scale_fill_viridis_c() +
  # scale_fill_gradientn(colours=viridis(7),na.value = "transparent",
  #                      breaks=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5),
  #                      limits=c(0,1))+
  # guides(alpha = FALSE, 
  #        size = FALSE) + 
  guides(alpha = FALSE) + 
  labs(
    x = "",
    y = "",
    title = "Desvio Padrão da Predição - Silte",
    size = "Silte Observado",
    fill = "Erro Padrão"
  ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

canit_dm <- gridExtra::grid.arrange(areia_std_err_dm,argila_std_err_dm,silte_std_err_dm,
                                    ncol =2)
