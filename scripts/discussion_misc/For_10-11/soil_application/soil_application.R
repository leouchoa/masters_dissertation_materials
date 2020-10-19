library(BivMaternEstim)
library(gstat)
library(sp)
library(ggplot2)
library(tidyr)
library(ggmap)

#krig_grid of leng 150 is good


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

initial_point <- c(1,1,0.5,0)
nus_vec <- c(1,1)
nug_vec <- c(0.5,0.5)

soil_dts_fit <- fit_biwm(
  obs_matrix = soil_dts@data,
  coords_matrix = coordinates(soil_dts),
  theta0 = initial_point,
  nus = nus_vec,
  nug_vec = nug_vec
  )

coords_range <- apply(coordinates(soil_dts),2,range)

krig_grid <-
  expand.grid(
    seq(coords_range[1,1] - 0.05,coords_range[2,1] + 0.09,length.out = 150),
    seq(coords_range[1,2] - 0.05,coords_range[2,2] + 0.09,length.out = 150)
  )

# krig_grid_full <- 
#   expand.grid(
#     seq(-21.1,-17.6,length.out = 50),
#     seq(-49.9,-42.9,length.out = 50)
#   )

soil_dts_preds <- compositional_biwm_krig(
  biwm_fit = soil_dts_fit,
  krig_locations = krig_grid,
  fit_locations = coordinates(soil_dts),
  obs_matrix = soil_dts@data,
  nus = nus_vec,
  nug_vec
  )


# soil_dts_preds_full <- compositional_biwm_krig(
#   biwm_fit = soil_dts_fit,
#   krig_locations = krig_grid_full,
#   fit_locations = coordinates(soil_dts),
#   obs_matrix = soil_dts@data,
#   nus = nus_vec,
#   nug_vec
# )

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

# soil_dts_preds_plot


# soil_dts_preds_full_plot <- soil_dts_preds_full %>% 
#   gather("key","Valor",1:3) %>% 
#   ggplot(aes(coord_x,coord_y, fill = Valor)) +
#   geom_tile() +
#   facet_wrap(~key, labeller = 
#                labeller(key =c(comp_01 = "Composição 1",
#                                comp_02 = "Composição 2",
#                                comp_03 = "Composição 3")),
#              ncol = 1) + 
#   labs(x = "", y = "") +
#   scale_fill_viridis_c()
# 
# soil_dts_preds_full_plot



# --------- Prediction Map --------------------

ctb_map_image <- readRDS("ctb_map_image_file")
ctb_map_image_unique_loc_constrained <- readRDS("ctb_map_image_unique_loc_constrained.rda")

ctb_map_image_unique_loc_constrained_v2 <- readRDS("ctb_map_image_unique_loc_constrained_v2.rda")


soil_dts_map_plot_v2 <- ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")), 
             mapping = aes(x = coord_x,y = coord_y))


names(soil_dts_preds) <- c("areia","argila","silte","coord_x","coord_y")

pred_map_areia <- soil_dts_map_plot_v2 + 
  geom_tile(data = soil_dts_preds[,c(1,4,5)], 
            aes(coord_x,coord_y, fill = areia, alpha = 0.01)
  ) + 
  scale_fill_viridis_c() +
  guides(alpha = FALSE) + 
  # guides(fill=guide_legend(title="Proporção \n de Areia")) +
  # guides(fill=guide_legend(title="")) +
  labs(
    x = "",
    y = "",
    title = "Predição Composicional de Areia"
  ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
        )


pred_map_argila <- soil_dts_map_plot_v2 + 
  geom_tile(data = soil_dts_preds[,c(2,4,5)], 
            aes(coord_x,coord_y, fill = argila, alpha = 0.1)
  ) + 
  scale_fill_viridis_c() +
  guides(alpha = FALSE) + 
  labs(
    x = "",
    y = "",
    title = "Predição Composicional de Argila"
  )+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

pred_map_silte <- soil_dts_map_plot_v2 + 
  geom_tile(data = soil_dts_preds[,c(3,4,5)], 
            aes(coord_x,coord_y, fill = silte, alpha = 0.2)
  ) + 
  scale_fill_viridis_c() +
  guides(alpha = FALSE) + 
  labs(
    x = "",
    y = "",
    title = "Predição Composicional de Silte"
  )+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )


b <- gridExtra::grid.arrange(pred_map_areia,pred_map_argila,pred_map_silte,nrow = 3)

ggsave("comp_pred_nrow_3_no_axis.pdf",plot = b,width = 8.27,height = 11.69)




# --------- Prediction Map with SIZE --------------------


pred_map_areia_with_size <- ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")),mapping = aes(x = coord_x,y = coord_y,size = alr_transformed_with_locations_UNIQUE$prop1_areia)) + 
  geom_tile(data = soil_dts_preds[,c(1,4,5)], 
            aes(coord_x,coord_y, fill = areia, alpha = 0.01)
  ) + 
  scale_fill_viridis_c() +
  guides(alpha = FALSE,size = FALSE) + 
  # guides(fill=guide_legend(title="Proporção \n de Areia")) +
  # guides(fill=guide_legend(title="")) +
  labs(
    x = "",
    y = "",
    title = "Predição Composicional de Areia"
  ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

ggsave("pred_map_areia_with_size.pdf",plot = pred_map_areia_with_size,width = 8.27,height = 11.69)

pred_map_argila_with_size <- ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")),mapping = aes(x = coord_x,y = coord_y,size = alr_transformed_with_locations_UNIQUE$prop2_argila)) + 
  geom_tile(data = soil_dts_preds[,c(2,4,5)], 
            aes(coord_x,coord_y, fill = argila, alpha = 0.01)
  ) + 
  scale_fill_viridis_c() +
  guides(alpha = FALSE,size = FALSE) + 
  # guides(fill=guide_legend(title="Proporção \n de Areia")) +
  # guides(fill=guide_legend(title="")) +
  labs(
    x = "",
    y = "",
    title = "Predição Composicional de Argila"
  ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

ggsave("pred_map_argila_with_size.pdf",plot = pred_map_argila_with_size,width = 8.27,height = 11.69)

pred_map_silte_with_size <- ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")),mapping = aes(x = coord_x,y = coord_y,size = alr_transformed_with_locations_UNIQUE$prop3_silte)) + 
  geom_tile(data = soil_dts_preds[,c(3,4,5)], 
            aes(coord_x,coord_y, fill = silte, alpha = 0.01)
  ) + 
  scale_fill_viridis_c() +
  guides(alpha = FALSE,size = FALSE) + 
  # guides(fill=guide_legend(title="Proporção \n de Areia")) +
  # guides(fill=guide_legend(title="")) +
  labs(
    x = "",
    y = "",
    title = "Predição Composicional de Silte"
  ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

ggsave("pred_map_silte_with_size.pdf",plot = pred_map_silte_with_size,width = 8.27,height = 11.69)


b_with_size <- gridExtra::grid.arrange(pred_map_areia_with_size,pred_map_argila_with_size,pred_map_silte_with_size,nrow = 3)

ggsave("comp_pred_nrow_3_no_axis_with_size.pdf",plot = b_with_size,width = 8.27,height = 11.69)



# --------- Prediction Map same scale with SIZE --------------------

pred_map_areia_same_scale_with_size <- 
  ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(
    data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")),
    mapping = aes(x = coord_x,
                  y = coord_y,
                  size = alr_transformed_with_locations_UNIQUE$prop1_areia
                  )
    ) + 
  geom_tile(data = soil_dts_preds[,c(1,4,5)], 
            aes(coord_x,coord_y, fill = areia, alpha = 0.2)
  ) + 
  scale_fill_gradientn(colours=viridis(7),na.value = "transparent",
                       breaks=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5),
                       limits=c(0,1))+
  guides(alpha = FALSE, size = FALSE) +   
  labs(
    x = "",
    y = "",
    title = "Predição Composicional de Areia",
    size = "Areia"
  ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )


pred_map_argila_same_scale_with_size <- 
  ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(
    data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")),
    mapping = aes(x = coord_x,
                  y = coord_y,
                  size = alr_transformed_with_locations_UNIQUE$prop2_argila
    )
  ) +
  geom_tile(data = soil_dts_preds[,c(2,4,5)], 
            aes(coord_x,coord_y, fill = argila, alpha = 0.2)
  ) + 
  # scale_fill_viridis_c(begin = 0,end = 1) +
  scale_fill_gradientn(colours=viridis(7),na.value = "transparent",
                       breaks=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5),
                       limits=c(0,1))+
  guides(alpha = FALSE, size = FALSE) + 
  labs(
    x = "",
    y = "",
    title = "Predição Composicional de Argila",
    size = "Argila"
  ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )


pred_map_silte_same_scale_with_size <- 
  ggmap(ctb_map_image_unique_loc_constrained_v2) + 
  geom_point(
    data = setNames(as.data.frame(soil_dts@coords),c("coord_x","coord_y")),
    mapping = aes(x = coord_x,
                  y = coord_y,
                  size = alr_transformed_with_locations_UNIQUE$prop3_silte
    )
  ) + 
  geom_tile(data = soil_dts_preds[,c(3,4,5)], 
            aes(coord_x,coord_y, fill = silte, alpha = 0.2)
  ) + 
  # scale_fill_viridis_c(begin = 0,end = 1) +
  scale_fill_gradientn(colours=viridis(7),na.value = "transparent",
                       breaks=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5),
                       limits=c(0,1))+
  guides(alpha = FALSE, size = FALSE) + 
  labs(
    x = "",
    y = "",
    title = "Predição Composicional de Silte",
    size = "Silte"
  ) +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

comp_pred_nrow_3_no_axis_same_scale_with_size <- gridExtra::grid.arrange(pred_map_areia_same_scale_with_size,pred_map_argila_same_scale_with_size,pred_map_silte_same_scale_with_size)

ggsave("comp_pred_nrow_3_no_axis_same_scale_with_size.pdf",plot = comp_pred_nrow_3_no_axis_same_scale_with_size,width = 8.27,height = 11.69)


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
