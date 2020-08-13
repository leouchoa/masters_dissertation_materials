create_sigma2_1_deriv <- function(sigmas,
                                  as,
                                  rho,
                                  coords_matrix,
                                  combined = TRUE){
  
  d <- dist(coords_matrix)
  a <- as[3] #ffs, what a shitty code i've made
  
  D_11_sigma2_2 <- matrix(0,nrow(coords_matrix),nrow(coords_matrix))
  
  
  D_21_sigma2_2 <- 0.5 * rho * sqrt(sigmas[1]/ sigmas[2]) * 
    BivMaternEstim:::matern_cov_wrapper(d,
                                        a = a,
                                        nu = nus[3]
                                        )
  
  
  D_22_sigma2_2 <- sigmas[1] * BivMaternEstim:::matern_cov_wrapper(d,
                                                                   a = a,
                                                                   nu = nus[2]
  )
  
  
  # Matrix::isSymmetric(Deriv_matrix_sigma2_2)
  
  if(combined){
    return(
      cbind(rbind(D_11_sigma2_2, D_21_sigma2_2), rbind(D_21_sigma2_2, D_22_sigma2_2))
    )
  }else return(out)
  
}



# ---- test it ----

# set.seed(123)
# coords_matrix <- matrix(runif(2*100), ncol = 2)
# d <- dist(coords_matrix)
# sigmas <- c(1,1)
# nus <- c(0.5,0.5);nus[3] <- mean(nus)
# as <- rep(2,3)
# rho <- 0.5
# 
# create_sigma2_1_deriv(sigmas,as,rho,coords_matrix)
# 
