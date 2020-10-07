soil_dts_fit_mock <- soil_dts_fit

# soil_dts_fit_mock$theta <- c(1,1,0.4,-0.1)
soil_dts_fit_mock$theta <- c(1,1,0.5,0)

initial_point <- c(1,1,0.5,0)

coords_range <- apply(coordinates(soil_dts),2,range)

krig_grid <- 
  expand.grid(
    seq(coords_range[1,1],coords_range[2,1],length.out = 30),
    seq(coords_range[1,2],coords_range[2,2],length.out = 30)
  )

soil_dts_preds <- compositional_biwm_krig(
  biwm_fit = soil_dts_fit_mock,
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

#thats how it should be, it makes much more sense

soil_dts_preds_plot + soil_dts_map_plot

# --------- New data points ----------


NEW_obs_soil_dts_preds_mock <- compositional_biwm_krig(
  biwm_fit = soil_dts_fit_mock,
  krig_locations = NEW_alr_transformed_with_locations_UNIQUE[,6:7],
  fit_locations = coordinates(soil_dts),
  obs_matrix = soil_dts@data,
  nus = nus_vec
)

# the error is large
colSums(NEW_obs_soil_dts_preds_mock[,1:3] - NEW_alr_transformed_with_locations_UNIQUE[,1:3])


# the map predictions makes much more sense, but the hist didnt improve a lot. I think it is misleading
apply(
  X = NEW_obs_soil_dts_preds_mock[,1:3] - NEW_alr_transformed_with_locations_UNIQUE[,1:3],
  MARGIN = 2,
  FUN = hist
)
