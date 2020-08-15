#' # TO-DO: Vector Field with Contour 
#' 
#' 1. Criar um grid para avaliar o gradiente do gráfico de contorno
#' 2. Avaliar o gradiente para o grid criado, adicionar o tamanho do passo
#' 3. Tomar cuidado com início e fim das flechas
#' 
#' ## Sources:
#' 
#' - https://dahtah.github.io/imager/gradient_field.html
#' - https://codereview.stackexchange.com/questions/93974/plotting-a-vector-field-in-r-ggplot
#' - https://stackoverflow.com/questions/14936504/vector-field-visualisation-r
#'
#' A função `eval_llike` avalia a log-verossimilhança para um vetor `theta_vec` fornecido.
#'
#' A função `plot_llike_contour` faz o gráfico de contorno para uma matriz em que cada observação é `theta_vec`.
#'
#' A função `construct_grid` é só pra deixar as coisas mais organizadas.
#'
#'   O código abaixo é importante para o funcionamento do codigo. Eu até posso consertar isso, mas acho que vai ficar pra depois.
#'
#'
#' set.seed(3)
#' n <- 40
#' coords <- matrix(runif(2*n), ncol = 2)
#' temp <- rnorm(2*n)
#'
#' nus_vec <- c(0.5,0.5)

eval_llike <- function(theta_vec,log_cd, coords){

  # S <- sigma_assembler_biwm(sigmas = theta_vec[1:2],
  #                           a = theta_vec[3],
  #                           rho = theta_vec[4],
  #                           nus = nus_vec,
  #                           coords_matrix = coords,
  #                           combined = TRUE)

  # log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)

  LLike_biwm(theta = theta_vec,
             nus = c(0.5,0.5),
             mu = colMeans(log_cd),
             coords_matrix = coords,
             obs_matrix = log_cd)

}

plot_llike_contour <- function(var_1,var_2,grid_df,obs, coords,true_param = NULL, n_bins = 30){


  eval_llike_results <- rep(0,nrow(grid_df))

  for(i in 1:nrow(grid_df)){

    eval_llike_results[i] <- eval_llike(
      as.matrix(grid_df[i,]),log_cd = obs,coords)

  }

  # ggplot(grid_df,
  #        aes_string(var_1, var_2, z = eval_llike_results)) + geom_contour()
  
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



plot_llike_contour_field <- function(var_1,var_2,grid_df,obs, coords,n_bins = 30,true_param = NULL,gamma_plot=0.005,...){
  
  
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

  
  arrows(x0 = grid_df[,var_1], y0 = grid_df[,var_2],
         x1 = grid_df[,var_1] + gamma_plot*grad_matrix[,var_1], y1 = grid_df[,var_2] + gamma_plot*grad_matrix[,var_2],
         length = 0.05, lwd = 2)

}


plot_llike_contour_field_fisher <- function(var_1,var_2,grid_df,obs, true_theta,coords,n_bins = 30,true_param = NULL,gamma_plot=0.005){
  
  source("fisher_info.R")
  
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
  
  
  fisher_info_matrix <- get_fisher_info(theta_vec = true_theta,
                                        nus = c(0.5,0.5),
                                        coords_matrix = coords)
  
  grad_matrix <- grad_matrix %*% solve(fisher_info_matrix)
  
  colnames(grad_matrix) <- c("sigma2_1","sigma2_2","a","rho")
  
  arrows(x0 = grid_df[,var_1], y0 = grid_df[,var_2],
         x1 = grid_df[,var_1] + gamma_plot*grad_matrix[,var_1], y1 = grid_df[,var_2] + gamma_plot*grad_matrix[,var_2],
         length = 0.05, lwd = 2)
 
  
   
}


# ---------- tests ------
# suppressPackageStartupMessages(library(BivMaternEstim))
# set.seed(123)
# 
# 
# 
# n <- 100
# coords <- matrix(runif(2*n), ncol = 2)
# temp <- rnorm(2*n)
# nus_vec <- c(0.5,0.5)
# 
# true_theta <- c(1,1,2,0.5)
# 
# S <- sigma_assembler_biwm(sigmas = true_theta[1:2], 
#                           a = true_theta[3], 
#                           rho = true_theta[4],
#                           nus = c(0.5, 0.5), 
#                           coords_matrix = coords, 
#                           combined = TRUE)
# 
# log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)
# 
# grid_02 <- construct_grid(
#   1,
#   1,
#   seq(0.5,3,length.out = 15),
#   seq(0.4,0.9,length.out = 15)
# )
# 
# par(mfrow = c(1,2))
# 
# plot_llike_contour_field_fisher("a","rho",grid_df = grid_02,obs = log_cd)
# plot_llike_contour_field("a","rho",grid_df = grid_02,obs = log_cd)
# 
# 
# plot_llike_contour_field_fisher("a","rho",grid_df = grid_02,obs = log_cd)
# 
# 
# grid_01 <- construct_grid(
#   seq(0.8,2.5,length.out = 15),
#   1,
#   seq(0.8,2.5,length.out = 15),
#   0.2)
# 
# plot_llike_contour_field("sigma2_1","a",grid_df = grid_01,obs = log_cd)
# 
# grid_03 <- construct_grid(
#   seq(0.25,1.5,length.out = 15),
#   seq(0.25,1.5,length.out = 15),
#   2,
#   0.5)
# 
# plot_llike_contour_field("sigma2_1","sigma2_2",grid_df = grid_03,obs = log_cd)
# 
# grid_04 <- construct_grid(
#   seq(0.25,1.5,length.out = 15),
#   1,
#   1,
#   seq(0.4,0.8,length.out = 15))
# 
# plot_llike_contour_field("sigma2_1","rho",grid_df = grid_04,obs = log_cd,gamma_plot = 0.0005)
