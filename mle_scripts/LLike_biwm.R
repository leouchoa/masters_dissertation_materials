# suppressPackageStartupMessages(library(fields))
source("../sigma_assembler_biwm.R")

#the order of params is: sigmas; as; rho; mu; nus

LLike_biwm <- function(theta,
                       nus,
                       mu,
                       coords_matrix,
                       obs_vec, 
                       parsimonious = TRUE,
                       ...){
  
  other_params <- list(...)
  
  if(parsimonious == TRUE){
    
    sigmas <- theta[1:2]
    as <- theta[3]
    rho <- theta[4]
    # mu <- theta[5]
    # nus <- theta[6:7]
    
    if(length(nus) != 2){
      stop("Smoothness parameter vector must be of length 2", 
           call. = FALSE)
    }
    
    if(length(unique(as)) != 1){
      stop("All shape paramaters must be the same", 
           call. = FALSE)
    }
    
    if(
      abs(rho) > ( sqrt(nus[1]*nus[2]) )/mean(nus)
    ){
      stop("Rho inserted does not define a valid covariance matrix",
           call. = FALSE)
    }
    
    #lazy (noob) preproc 
    nus[3] <- mean(nus[1:2])
    if(length(as) == 1){
      as <- rep(as,3)
    }
    
    autocov_matrix <- sigma_assembler_biwm(sigmas = sigmas,
                                           as = as,
                                           rho = rho,
                                           nus = nus,
                                           coords_matrix = coords_matrix)
    
  }
  
  else{
    
    sigmas <- theta[1:2]
    as <- theta[3:5]
    rho <- theta[6]
    # mu <- theta[7]
    # nus <- theta[8:10]
    
    autocov_matrix <- sigma_assembler_biwm(sigmas = sigmas,
                                           as = as,
                                           rho = rho,
                                           nus = nus,
                                           coords_matrix = coords_matrix)
  }
  
  # print(sigmas,as,rho)
  
  as.double(
    -1/2 *
      (
        log( det(autocov_matrix) ) +
          t( (obs_vec - mu) ) %*%
          solve(autocov_matrix, obs_vec - mu) +
          length(obs_vec) * log(2*pi)
      )
  )
  
}
