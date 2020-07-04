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

eval_llike <- function(theta_vec,log_cd){

  # S <- sigma_assembler_biwm(sigmas = theta_vec[1:2],
  #                           a = theta_vec[3],
  #                           rho = theta_vec[4],
  #                           nus = nus_vec,
  #                           coords_matrix = coords,
  #                           combined = TRUE)

  # log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)

  LLike_biwm(theta = theta_vec,
             nus = nus_vec,
             mu = colMeans(log_cd),
             coords_matrix = coords,
             obs_matrix = log_cd)

}

plot_llike_contour <- function(var_1,var_2,grid_df,obs,true_param = NULL, n_bins = 30){


  eval_llike_results <- rep(0,nrow(grid_df))

  for(i in 1:nrow(grid_df)){

    eval_llike_results[i] <- eval_llike(
      as.matrix(grid_df[i,]),log_cd = obs)

  }

  # ggplot(grid_df,
  #        aes_string(var_1, var_2, z = eval_llike_results)) + geom_contour()
  
  contour(
    unique(grid_df[,var_1]),
    unique(grid_df[,var_2]),
    matrix(eval_llike_results,ncol = 15),
    nlevels = n_bins,
    xlab = var_1,
    ylab = var_2)
    #main = "Contornos da Log-Verossimilhança")
  
  if( !is.null(true_param) ){
    points(true_param[1],true_param[2],pch = 8)
  }

}


plot_llike_contour_ggplot <- function(var_1,var_2,grid_df,obs,bins = NULL){
  
  
  eval_llike_results <- rep(0,nrow(grid_df))
  
  for(i in 1:nrow(grid_df)){
    
    eval_llike_results[i] <- eval_llike(
      as.matrix(grid_df[i,]),log_cd = obs)
    
  }
  
  ggplot(grid_df,
         aes_string(var_1, var_2, 
                    z = eval_llike_results)
         ) + 
    geom_contour(colour='black',bins = bins) +
    geom_dl(
      aes(label=..level..),
      method="bottom.pieces", 
      stat="contour"
    )
  
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
