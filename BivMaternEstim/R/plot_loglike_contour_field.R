#' Likelihood Contour Plot for the Bivariate (Wittle-)Matern Model.
#'
#' @export
#'
#' @examples
#'
#' n <- 70
#' coords <- matrix(runif(2*n), ncol = 2)
#' temp <- rnorm(2*n)
#' nus_vec <- c(0.5,0.5)
#' sample.int(123123,1)
#'
#'
#'
#' true_theta = c(1,1,2,0.5)
#' sigma2_1_grid <- seq(0.5,2,length.out = 15)
#' sigma2_2_grid <- seq(0.5,2,length.out = 15)
#' a_grid <- seq(0.5,3,length.out = 15)
#' rho_grid <- seq(0.25,0.75,length.out = 15)
#'
#' S <- sigma_assembler_biwm(sigmas = true_theta[1:2],
#'                           a = true_theta[3],
#'                           rho = true_theta[4],
#'                           nus_vector = nus_vec,
#'                           coords_matrix = coords,
#'                           combined = TRUE)
#'
#' log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)
#'
#' grid_01 <- construct_grid(
#'   sigma2_1_grid,
#'   1,
#'   a_grid,
#'   0.2)
#'
#' grid_02 <- construct_grid(
#'   1,
#'   1,
#'   a_grid,
#'   rho_grid
#' )
#'
#'
#' grid_03 <- construct_grid(
#'   sigma2_1_grid,
#'   sigma2_2_grid,
#'   2,
#'   0.5)
#'
#' grid_04 <- construct_grid(
#'   sigma2_1_grid,
#'   1,
#'   1,
#'   rho_grid
#' )
#'
#'
#' plot_loglike_contour_ggplot("sigma2_1","a",grid_df = grid_01,obs = log_cd,coords = coords,nbins = 30,gamma_plot =1e3,true_param = true_theta[c(1,3)])
#'
#'
#' plot_loglike_contour_ggplot("a","rho",grid_df = grid_02,obs = log_cd,coords = coords,nbins = 30,gamma_plot =1e3,true_param = true_theta[c(3,4)])
#'
#'
#'
#' plot_loglike_contour_ggplot("sigma2_1","sigma2_2",grid_df = grid_03,obs = log_cd,coords = coords,nbins = 30,gamma_plot =1e3,true_param = true_theta[c(1,2)])
#'
#'
#'
#' plot_loglike_contour_ggplot("sigma2_1","rho",grid_df = grid_04,obs = log_cd,coords = coords,nbins = 30,gamma_plot =1e3,true_param = true_theta[c(1,4)])
#'
#'


plot_loglike_contour_ggplot <- function(var_1,var_2,grid_df,obs, coords,nus_vector,true_param = NULL,nbins = 30,gamma_plot = 1){


  eval_loglike_results <- rep(0,nrow(grid_df))

  for(i in 1:nrow(grid_df)){

    eval_loglike_results[i] <- eval_loglike(
      as.matrix(grid_df[i,]),log_cd = obs,coords)

  }

  grad_matrix <-
    matrix(nrow = nrow(grid_df),ncol = ncol(grid_df))


  for(i in seq(nrow(grid_df))){

    grad_matrix[i,] <-  block_LLike_biwm_grad(
      theta = as.matrix(grid_df[i,]),
      nus = nus_vector,
      coords_matrix = coords,
      obs_matrix = obs,
      mu = colMeans(obs)
    )
  }

  colnames(grad_matrix) <- c("sigma2_1","sigma2_2","a","rho")

  contour_plot <- ggplot(grid_df,
         aes_string(var_1, var_2, z = eval_loglike_results)
  ) +
    geom_contour(bins = nbins) +
    geom_segment(
      aes(
        x = grid_df[,var_1],
        y = grid_df[,var_2],
        xend = grid_df[,var_1] + gamma_plot*grad_matrix[,var_1],
        yend = grid_df[,var_2] + gamma_plot*grad_matrix[,var_2]
      ),
      arrow = arrow(length = unit(0.15,"cm"))
    ) + theme_minimal()

  if(!is.null(true_param)){
    return(
      contour_plot +
        geom_point(
          aes(x = true_param[1],y = true_param[2], color = "red"),
          shape = 8,
          size = 3
        )
    )
  }else{
    return(
      contour_plot
    )
  }



}

