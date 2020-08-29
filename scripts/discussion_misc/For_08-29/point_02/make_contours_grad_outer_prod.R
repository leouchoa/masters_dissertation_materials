make_contours_grad_outer_prod <- function(var_1,var_2,grid_df,obs, true_theta,coords,n_bins = 30,true_param = NULL,gamma_plot=0.005){
  
  
  eval_llike_results <- rep(0,nrow(grid_df))
  
  for(i in 1:nrow(grid_df)){
    
    eval_llike_results[i] <- eval_llike(
      as.matrix(grid_df[i,]),log_cd = obs,coords)
    
  }
  
  contour(
    unique(grid_df[,var_1]),
    unique(grid_df[,var_2]),
    matrix(eval_llike_results,ncol = 15),
    nlevels = n_bins,
    xlab = var_1,
    ylab = var_2,
    main = "Contornos da Log-Verossimilhança")
  
  if( !is.null(true_param) ){
    points(true_param[1],true_param[2],pch = 8)
  }
  
  grad_matrix <- 
    matrix(nrow = nrow(grid_df),ncol = ncol(grid_df))
  
  
  for(i in seq(nrow(grid_df))){
    
    grad_matrix[i,] <-  BivMaternEstim:::block_LLike_biwm_grad(
      theta = as.matrix(grid_df[i,]),
      nus = c(0.5, 0.5),
      coords_matrix = coords,
      obs_matrix = obs,
      mu = colMeans(obs)
    )
  }
  
  colnames(grad_matrix) <- c("sigma2_1","sigma2_2","a","rho")
  
  
  

  
  grad_matrix <- grad_matrix %*% solve(fisher_info_matrix)
  
  colnames(grad_matrix) <- c("sigma2_1","sigma2_2","a","rho")
  
  arrows(x0 = grid_df[,var_1], y0 = grid_df[,var_2],
         x1 = grid_df[,var_1] + gamma_plot*grad_matrix[,var_1], y1 = grid_df[,var_2] + gamma_plot*grad_matrix[,var_2],
         length = 0.05, lwd = 2)
  
  
  
}