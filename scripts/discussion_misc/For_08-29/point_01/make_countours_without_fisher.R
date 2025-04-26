#' It's necessary to the script `contourField_funs.R` to be in the same folder as the file contains complementary functions
#' 
#' 
#' 
#' IT DOESNT WORK IF THE GRADIENT IS ALREADY PRE-MULTIPLIED BY `gamma_plot_val` because it will be really small. A **not** possible workaround its to pass `gamma_plot_val` with inverse scale 

make_countours_without_fisher <- function(true_theta = 
                                          c(1,1,2,0.5), 
                                          n =70,
                                          gamma_plot_val = 0.0005
                                          ){

  suppressPackageStartupMessages(library(BivMaternEstim))
  set.seed(123)
  
  source("contourField_funs.R")



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

    
    
    
    
    
  plot_llike_contour_field("sigma2_1","a",grid_df = grid_01,obs = log_cd,coords = coords,gamma_plot = gamma_plot_val)
  
  
  grid_02 <- construct_grid(
    1,
    1,
    seq(1,3,length.out = 15),
    seq(0.25,0.75,length.out = 15)
  )
  
  plot_llike_contour_field("a","rho",grid_df = grid_02,obs = log_cd,coords = coords, gamma_plot = gamma_plot_val)

  
  
  
  
  
  grid_03 <- construct_grid(
    seq(0.25,1.5,length.out = 15),
    seq(0.25,1.5,length.out = 15),
    2,
    0.5)


  plot_llike_contour_field("sigma2_1","sigma2_2",grid_df = grid_03,obs = log_cd,gamma_plot = gamma_plot_val)

  
  
  
  grid_04 <- construct_grid(
    seq(0.25,2,length.out = 15),
    1,
    1,
    seq(0.25,0.75,length.out = 15)
  )

  plot_llike_contour_field("sigma2_1","rho",grid_df = grid_04,obs = log_cd,coords = coords,gamma_plot = gamma_plot_val)


}


# ----- Test -----

make_countours_without_fisher()

gamma_vec <- c(0.005,0.0005,0.00005,0.000005,0.0000005)

for(i in gamma_vec[1]){
  
  file_name = paste0(
    "contours_gamma_plot_",
    gamma_vec[i],
    ".pdf"
    )
  
  pdf("rplot.pdf", width = 8.27, height = 11.69)
  
  make_countours_without_fisher(gamma_plot_val = gamma_vec[i])
  
  dev.off()
}


