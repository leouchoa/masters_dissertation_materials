#' for now this `soil_dts_fit` is needed from soil_application.R

# Step 01: krigging
alr_dts_preds <- 
  setNames(
    compositional_biwm_krig(
      biwm_fit = soil_dts_fit,
      krig_locations = krig_grid,
      fit_locations = coordinates(soil_dts),
      obs_matrix = soil_dts@data,
      nus = nus_vec,
      nug_vec,
      return_comp = FALSE
    ),
    c("h_1","h_2","coord_x","coord_y")
  )


# Step 02: transform krig estimates

head(soil_dts_preds)

step_02_df <- 
  cbind(
    data.frame(
      hz_1 = exp(alr_dts_preds$h_1),
      hz_2 = exp(alr_dts_preds$h_2),
      hz_3 = 1/( 1 + exp(alr_dts_preds$h_1) + exp(alr_dts_preds$h_2) )
    ),
    alr_dts_preds[,c(3:4)]
  )



# Step 03: cond var

sigma_bar <- BivMaternEstim:::conditional_cov_mat(
  loc_new = krig_grid,
  loc_obs = soil_dts@coords,
  sigmas = soil_dts_fit$theta[1:2],
  a = soil_dts_fit$theta[3],
  soil_dts_fit$theta[4],
  nus = nus_vec,
  nug_vec = nug_vec
)





# Step 04: gradient


extract_block_sigbar <- function(sig_bar){
  
  limits <- dim(sig_bar)/2
  
  block_sig11 <- sig_bar[1:limits[1],1:limits[2]]
  block_sig12 <- sig_bar[1:limits[1],(limits[2]+1):ncol(sig_bar)]
  block_sig21 <- sig_bar[(limits[1]+1):nrow(sig_bar),1:limits[2]]
  block_sig22 <- sig_bar[(limits[1]+1):nrow(sig_bar),(limits[2]+1):ncol(sig_bar)]
  
  return(
    list(
      block_sig11 = block_sig11,
      block_sig12 = block_sig12,
      block_sig21 = block_sig21,
      block_sig22 = block_sig22
    )
  )
}



extract_individual_sigbar <- function(iter,extract_block_sigbar_res){
  
  aux_mat <- matrix(0,2,2)
  
  aux_mat[1,1] <- extract_block_sigbar_res$block_sig11[iter,iter]
  aux_mat[1,2] <- extract_block_sigbar_res$block_sig12[iter,iter]
  aux_mat[2,1] <- extract_block_sigbar_res$block_sig21[iter,iter]
  aux_mat[2,2] <- extract_block_sigbar_res$block_sig22[iter,iter]
  
  return(
    aux_mat
  )
}


create_delta_grad_mat <- function(mu){
  
  #ffs
  if(is.data.frame(mu)){
    mu <- c(mu[[1]],mu[[2]])
  }
  
  dm_grad_mat <- matrix(0,2,3)
  
  common_denominator <- (1 + exp(mu[1]) + exp(mu[2]))^2
  
  dm_grad_mat[1,1] <- exp(mu[1]) * (1 + exp(mu[2]))
  dm_grad_mat[1,2] <- - exp(mu[1]) * exp(mu[2])
  dm_grad_mat[1,3] <- - exp(mu[1]) 
  
  dm_grad_mat[2,1] <- - exp(mu[1]) * exp(mu[2])
  dm_grad_mat[2,2] <- exp(mu[2]) * (1 + exp(mu[1]))
  dm_grad_mat[2,3] <- - exp(mu[2]) 
  
  
  
  return(
    dm_grad_mat/common_denominator
  )
}


#NOW I NEED TO ITERATE OVER LENGTH(DIAG(SIGMA_BAR)) AND MAKE THE PROCEDURE DESCRIBE IN THE DISSERTATION TO GET IT

# Step 05: calculate approx cov mat

get_approx_alrInv_cond_covvar <- function(sigma_bar,
                                          krig_estimates_mat,
                                          coords_matrix){
  
  N <- nrow(sigma_bar) / 2
  
  sig_bar_blocks <- extract_block_sigbar(sig_bar = sigma_bar)
  
  approx_alrInv_cond_covvar <- matrix(0,N,3)
  
  for(i in 1:N){
    single_sig_bar <- extract_individual_sigbar(i,sig_bar_blocks)
    
    approx_mean <- create_delta_grad_mat(krig_estimates_mat[i,])
    
    approx_alrInv_cond_covvar[i,] <- 
      diag(t(approx_mean) %*% single_sig_bar %*% approx_mean)
  }
  
  return(
    cbind(
      setNames(
        as.data.frame(approx_alrInv_cond_covvar),
        c("comp_01_sd","comp_02_sd","comp_03_sd")),
      coords_matrix
    )
  ) 
  
}


# debugonce(get_approx_alrInv_cond_covvar)

sigbar_hat_approx <- get_approx_alrInv_cond_covvar(sigma_bar = sigma_bar,
                              krig_estimates_mat = alr_dts_preds[,1:2],
                              coords_matrix = soil_dts_preds[,4:5]
                              )


sigbar_hat_approx %>% 
  gather("key","Valor",1:3) %>% 
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller =
               labeller(key =c(comp_01_sd = "Composição 1",
                               comp_02_sd = "Composição 2",
                               comp_03_sd = "Composição 3"))) +
  labs(x = "", y = "") +
  scale_fill_viridis_c() +
  labs(
    title = "Erro Padrão da Predição Composicional"
    # subtitle = paste("Pontos Amostrados:", 80)
  )


ggmap(ctb_map_image_unique_loc_constrained_v2) + 
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

ggmap(ctb_map_image_unique_loc_constrained_v2) + 
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

ggmap(ctb_map_image_unique_loc_constrained_v2) + 
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
