library(ggmap)
alr_transformed_with_locations_UNIQUE <- read.csv("alr_transformed_with_locations_UNIQUE.csv", row.names=NULL)


coords_range <- apply(alr_transformed_with_locations_UNIQUE[,c("coord_x","coord_y")],2,range)

# krig_grid <-
#   expand.grid(
#     seq(coords_range[1,1],coords_range[2,1],length.out = 50),
#     seq(coords_range[1,2],coords_range[2,2],length.out = 50)
#   )

height_unique <- unname(coords_range[2,2] - coords_range[1,2])
width_unique <- unname(coords_range[2,1] - coords_range[1,1])

sac_borders_unique_loc_constrained <- 
  c(bottom  = unname(coords_range)[1,2]  - 0.1 * height_unique, 
    top     = unname(coords_range)[2,2]  + 0.1 * height_unique,
    left    = unname(coords_range)[1,1] - 0.1 * width_unique,
    right   = unname(coords_range)[2,1] + 0.1 * width_unique)

sac_borders_unique_loc_constrained_v2 <- 
  c(bottom  = unname(coords_range)[1,2]  - 0.01 * height_unique, 
    top     = unname(coords_range)[2,2]  + 0.01 * height_unique,
    left    = unname(coords_range)[1,1] - 0.01 * width_unique,
    right   = unname(coords_range)[2,1] + 0.01 * width_unique)


ctb_map_image_unique_loc_constrained <- get_stamenmap(sac_borders_unique_loc_constrained, zoom = 7, maptype = "terrain")

ctb_map_image_unique_loc_constrained_v2 <- get_stamenmap(sac_borders_unique_loc_constrained_v2, zoom = 7, maptype = "terrain")


saveRDS(ctb_map_image_unique_loc_constrained,
        "ctb_map_image_unique_loc_constrained.rda")

saveRDS(ctb_map_image_unique_loc_constrained_v2,
        "ctb_map_image_unique_loc_constrained_v2.rda")

