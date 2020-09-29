gen_initial_points <- function(true_params){
  
  sigmas_initial_point <- rchisq(2,true_params[1:2])
  a_initial_point <- rchisq(1,true_params[3])
  rho_initial_point <- rnorm(1,true_params[4],0.3)
  rho_initial_point <- ifelse(rho_initial_point <= -1,-0.99,rho_initial_point)
  rho_initial_point <- ifelse(rho_initial_point >= 1,0.99,rho_initial_point)
  return(
    c(
      sigmas_initial_point,
      a_initial_point,
      rho_initial_point
    )
  )
}