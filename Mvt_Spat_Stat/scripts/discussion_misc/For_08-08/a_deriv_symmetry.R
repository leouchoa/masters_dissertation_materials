#' I make `M_dash_nu_k` as a function argument because in the future we're gonna call it from `block_LLike_biwm_grad` and so we already have `M_dash_nu_k` computed

create_a_deriv <- function(M_dash_nu_1,
                           M_dash_nu_2,
                           M_dash_nu_3,
                           sigmas,
                           as,
                           rho,
                           combined = TRUE){
  
  
  D_11_a <- sigmas[1] * M_dash_nu_1
  
  D_21_a <- rho * sqrt(sigmas[2]/ sigmas[1]) * M_dash_nu_3
  
  D_22_a <- sigmas[2] * M_dash_nu_2
  
  out <- cbind(rbind(D_11_a, D_21_a), rbind(t(D_21_a), D_22_a))
  
  # range(Deriv_matrix_a - t(Deriv_matrix_a))
  
  if(combined){
    return(
      cbind(rbind(D_11_a, D_21_a), rbind(t(D_21_a), D_22_a))
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
# M_dash_nu_1 <- BivMaternEstim:::matern_deriv(a = as[1], 
#                                              nu = nus[1], 
#                                              coords_matrix = coords_matrix)
# 
# M_dash_nu_2 <- BivMaternEstim:::matern_deriv(a = as[2], 
#                                              nu = nus[2], 
#                                              coords_matrix = coords_matrix)
# 
# M_dash_nu_3 <- BivMaternEstim:::matern_deriv(a = as[3], 
#                                              nu = nus[3], 
#                                              coords_matrix = coords_matrix)
# 
# create_a_deriv(M_dash_nu_1,M_dash_nu_2,M_dash_nu_3,sigmas,as,nus)
# 
