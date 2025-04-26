nus_vec <- c(1,1)


sigma2_1_grid <- seq(1,4,length.out = 15)
sigma2_2_grid <- seq(0.5,2,length.out = 15)
a_grid <- seq(46,56,length.out = 15)
rho_grid <- seq(-0.25,0.4,length.out = 15)

S <- sigma_assembler_biwm(sigmas = soil_dts_fit$theta[1:2],
                          a = soil_dts_fit$theta[3],
                          rho = soil_dts_fit$theta[4],
                          nus = nus_vec,
                          coords_matrix = coordinates(soil_dts),
                          combined = TRUE)


grid_01 <- construct_grid(
  sigma2_1_grid,
  1,
  a_grid,
  0.2)

plot_loglike_contour_ggplot("sigma2_1","a",grid_df = grid_01,obs = soil_dts@data,coords = coordinates(soil_dts),nbins = 30,gamma_plot =1e3,nus_vector = nus_vec)
