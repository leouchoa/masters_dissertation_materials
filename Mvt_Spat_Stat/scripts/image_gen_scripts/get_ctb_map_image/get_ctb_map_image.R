library(readxl)
library(ggmap)

ctb_observacao <- 
  read_excel(
    "ctb0562/ctb0562-observacao.xlsx", 
    col_types = c("text", "numeric", "text", 
                  "date", "skip", "numeric", "numeric", 
                  "skip", "skip", "skip", "skip", "text", 
                  "text", "text", "text", "text", "text", 
                  "skip", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "numeric", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text"))[-(1:2),]


coord_y_max <- max(ctb_observacao$coord_y,na.rm = TRUE)
coord_y_min <- min(ctb_observacao$coord_y,na.rm = TRUE)
coord_x_max <- max(ctb_observacao$coord_x,na.rm = TRUE)
coord_x_min <- min(ctb_observacao$coord_x,na.rm = TRUE)

height <- coord_y_max - coord_y_min
width <- coord_x_max - coord_x_min
sac_borders <- c(bottom  = coord_y_min  - 0.1 * height, 
                 top     = coord_y_max  + 0.1 * height,
                 left    = coord_x_min - 0.1 * width,
                 right   = coord_x_max + 0.1 * width)

ctb_map_image <- get_stamenmap(sac_borders, zoom = 7, maptype = "terrain")

saveRDS(ctb_map_image,"ctb_map_image_file")

