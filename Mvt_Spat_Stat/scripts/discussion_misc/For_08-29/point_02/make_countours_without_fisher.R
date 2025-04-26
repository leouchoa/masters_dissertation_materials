make_countours_without_fisher <- function(true_theta = 
                                          c(1,1,2,0.5), 
                                          n =70
                                          ){

  suppressPackageStartupMessages(library(BivMaternEstim))
  set.seed(123)



  coords <- matrix(runif(2*n), ncol = 2)
  temp <- rnorm(2*n)
  nus_vec <- c(0.5,0.5)



  S <- sigma_assembler_biwm(sigmas = true_theta[1:2],
                            a = true_theta[3],
                            rho = true_theta[4],
                            nus = c(0.5, 0.5),
                            coords_matrix = coords,
                            combined = TRUE)

  log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)

  par(mfrow = c(2,2), pty = "s")
  

    grid_01 <- construct_grid(
    seq(0.5,2,length.out = 15),
    1,
    seq(0.5,3,length.out = 15),
    0.2)

    
    
    
    
    
  plot_llike_contour_field("sigma2_1","a",grid_df = grid_01,obs = log_cd,coords = coords)
  
  
  grid_02 <- construct_grid(
    1,
    1,
    seq(1,3,length.out = 15),
    seq(0.25,0.75,length.out = 15)
  )
  
  plot_llike_contour_field("a","rho",grid_df = grid_02,obs = log_cd,coords = coords, gamma_plot = 0.0005)

  
  
  
  
  
  grid_03 <- construct_grid(
    seq(0.25,1.5,length.out = 15),
    seq(0.25,1.5,length.out = 15),
    2,
    0.5)


  plot_llike_contour_field("sigma2_1","sigma2_2",grid_df = grid_03,obs = log_cd,coords = coords)

  
  
  
  grid_04 <- construct_grid(
    seq(0.25,2,length.out = 15),
    1,
    1,
    seq(0.25,0.75,length.out = 15)
  )

  plot_llike_contour_field("sigma2_1","rho",grid_df = grid_04,obs = log_cd,coords = coords,gamma_plot = 0.0005)


}


# ----- Test -----

make_countours_without_fisher()

