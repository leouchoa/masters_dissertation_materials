suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_minimal())

eval_llike <- function(theta_vec,log_cd){
  
  LLike_biwm(theta = theta_vec,
             nus = nus_vec,
             mu = colMeans(log_cd),
             coords_matrix = coords,
             obs_matrix = log_cd)
  
}


plot_univariate_llike <- function(param,grid_df,obs,true_param = NULL,my_title = NULL){
  
  
  eval_llike_results <- rep(0,nrow(grid_df))
  
  for(i in 1:nrow(grid_df)){
    
    eval_llike_results[i] <- eval_llike(
      as.matrix(grid_df[i,]),log_cd = obs)
    
  }
  
  aux_df <- setNames(
    data.frame(grid_df[,param],eval_llike_results),c(param,"log_likelihood")
  )
  
  ggplot(aux_df,aes_string(param,"log_likelihood")) +
    geom_line() +
    geom_vline(xintercept=true_param) +
    labs(title = my_title)
  
  
}

construct_grid <- function(sigma2_1_vec,
                           sigma2_2_vec,
                           a_vec,
                           rho_vec){
  
  setNames(
    expand.grid(sigma2_1_vec,
                sigma2_2_vec,
                a_vec,
                rho_vec),
    c("sigma2_1","sigma2_2","a","rho"))
  
}


plot_both_zoom_llike <- function(not_zoomed,zoomed,param,n_col = NULL){
  gridExtra::grid.arrange(
    not_zoomed + 
      labs(title = paste("Not zoomed", param, "Log-like")),
    zoomed + 
      labs(title = paste("Zoomed", param, "Log-like")),
    ncol = n_col
  )
}

