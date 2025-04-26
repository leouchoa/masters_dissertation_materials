create_rho_deriv <- function(sigmas,
                             as,
                             rho,
                             coords_matrix,
                             combined = TRUE){
  
  
  a <- as[3] #ffs, what a shitty code i've made
  
  
  D_11_rho <- matrix(0,nrow(coords_matrix),nrow(coords_matrix))
  
  #i was receiving some strange errors because i was pulling d from another environment, so i hardcoded it again
  d <- dist(coords_matrix)
  
  D_21_rho <- sqrt(sigmas[2]/ sigmas[1]) * 
    BivMaternEstim:::matern_cov_wrapper(d,
                                        a = a,
                                        nu = nus[3]
    )
  
  D_22_rho <- matrix(0,nrow(coords_matrix),nrow(coords_matrix))
  
  
  # range(Deriv_matrix_rho - t(Deriv_matrix_rho))
  
  if(combined){
    return(
      cbind(rbind(D_11_rho, D_21_rho), rbind(t(D_21_rho), D_22_rho))
    )
  }else return(out)
  
  
}

# ---- test it ----
# 
# set.seed(123)
# coords_matrix <- matrix(runif(2*100), ncol = 2)
# d <- dist(coords_matrix)
# sigmas <- c(1,1)
# nus <- c(0.5,0.5);nus[3] <- mean(nus)
# as <- rep(2,3)
# rho <- 0.5
# 
# create_rho_deriv(sigmas,as,rho,coords_matrix)


