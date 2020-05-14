plot_rf_biwm <- function(rf_obj_as_df,x_grid,y_grid){
  
  library(ggplot2)
  
  my_grid <- expand.grid(x_grid,y_grid)
  
  gp_df <- data.frame(coord_x = my_grid$Var1,coord_y = my_grid$Var2, process_1 = rf_obj_as_df$variable1, process_2 = rf_obj_as_df$variable2)
  
  
  dts <- tidyr::gather(gp_df, "key","Valor",3:4)
  
  ggplot(dts, aes(coord_x,coord_y, fill = Valor)) +
    geom_tile() +
    facet_wrap(~key, labeller = 
                 labeller(key =c(process_1 = "Processo 1",
                                 process_2 = "Processo 2"))) + 
    labs(x = "", y = "") +
    scale_fill_viridis_c()
}