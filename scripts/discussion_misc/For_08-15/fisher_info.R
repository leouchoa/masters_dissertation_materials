get_fisher_info <- function(theta_vec,nus,coords_matrix){
  
  tr <- function(x){
    sum(diag(x))
    }
  
  # ------ Params Extraction -------  
  
  sigmas <- theta_vec[1:2]
  as <- rep(theta_vec[3],3)
  rho <- theta_vec[4]
  nus[3] <- mean(nus)
  
  d <- dist(coords_matrix)
  
  
  # ------ Inverse Calculation -------  
  
  sigma_matrix <- BivMaternEstim:::sigma_assembler_biwm(sigmas,
                                                        as[3],
                                                        rho,
                                                        nus,
                                                        coords_matrix,
                                                        combined = FALSE #must be FALSE
  )
  
  sigma_matrix_inv <- BivMaternEstim:::block_inv(sigma_matrix,combined = TRUE)
  
  # -------------- Derivative matrices ---------------
  
  
  #' # Pregame
  #' 
  #' 
  #' ## Creating Matern Derivative Matrices
  #' 
  #' 
  
  M_dash_nu_1 <- BivMaternEstim:::matern_deriv(a = as[1], 
                                               nu = nus[1], 
                                               coords_matrix = coords_matrix)
  
  M_dash_nu_2 <- BivMaternEstim:::matern_deriv(a = as[2], 
                                               nu = nus[2], 
                                               coords_matrix = coords_matrix)
  
  M_dash_nu_3 <- BivMaternEstim:::matern_deriv(a = as[3], 
                                               nu = nus[3], 
                                               coords_matrix = coords_matrix)
  
  #' # Caso em que $\theta = \sigma^2_1$ 
  #' 
  #' $B = \Sigma^{-1} \frac{\partial \Sigma}{\partial \sigma^2_1}$
  #' 
  #' 
  #' 
  #' 
  
  sigma2_1_deriv <- create_sigma2_1_deriv(sigmas,as,rho,coords_matrix)
  
  B_sigma2_1 <- tr(sigma_matrix_inv %*% sigma2_1_deriv)

  
  #' # Caso em que $\theta = \sigma^2_1$ 
  #' 
  #' $B = \Sigma^{-1} \frac{\partial \Sigma}{\partial \sigma^2_1}$
  #' 
  #' 
  #' 
  #' 
  
  sigma2_2_deriv <- create_sigma2_2_deriv(sigmas,as,rho,coords_matrix)
  
  B_sigma2_2 <- tr(sigma_matrix_inv %*% sigma2_2_deriv)
  
  
  #' # Caso em que $\theta = \rho$ 
  #' 
  #' $B = \Sigma^{-1} \frac{\partial \Sigma}{\partial \rho}$
  #' 
  #' 
  #' 
  #' 
  
  rho_deriv <- create_rho_deriv(sigmas,as,rho,coords_matrix)
  
  B_rho <- tr(sigma_matrix_inv %*% rho_deriv)
  
  
  
  #' # Caso em que $\theta = a$ 
  #' 
  #' $B = \Sigma^{-1} \frac{\partial \Sigma}{\partial \rho}$
  #' 
  #' 
  #' 
  #' 
  
  a_deriv <- create_a_deriv(M_dash_nu_1,M_dash_nu_2,M_dash_nu_3,sigmas,as,rho,coords_matrix)
  
  B_a <- tr(sigma_matrix_inv %*% a_deriv)
  
  
  # -------------- Fisher Information ---------------
  
  #' # Fisher information first line
  #' 
  #' 
  #' 
  
  #' ## Fisher information for $\sigma^2_1$
  
  F_11 <- 0.25 *(
    tr(sigma_matrix_inv %*% sigma2_1_deriv %*% 
         sigma_matrix_inv %*% sigma2_1_deriv) - 
      B_sigma2_1^2
      
  )
  
  #' ## Fisher information for $\sigma^2_1$ and $\sigma^2_2$
  
  F_12 <- 0.25 *(
    tr(sigma_matrix_inv %*% sigma2_1_deriv %*% 
         sigma_matrix_inv %*% sigma2_2_deriv) - 
      B_sigma2_1 * B_sigma2_2
    
  )
  
  #' ## Fisher information for $\sigma^2_1$ and $\rho$
  
  
  F_13 <- 0.25 *(
    tr(sigma_matrix_inv %*% sigma2_1_deriv %*% 
         sigma_matrix_inv %*% a_deriv) - 
      B_sigma2_1 * B_a
    
  )
  
  
  #' ## Fisher information for $\sigma^2_1$ and $a$
  
  
  F_14 <- 0.25 *(
    tr(sigma_matrix_inv %*% sigma2_1_deriv %*% 
         sigma_matrix_inv %*% rho_deriv) - 
      B_sigma2_1 * B_rho
    
  )
  
  #' # Fisher information second line
  #' 
  #' 
  #' ## Fisher information for $\sigma^2_2$
  
  F_22 <- 0.25 *(
    tr(sigma_matrix_inv %*% sigma2_2_deriv %*% 
         sigma_matrix_inv %*% sigma2_2_deriv) - 
      B_sigma2_2^2
    
  )
  
  
  #' ## Fisher information for $\sigma^2_2$ and $\rho$
  #' 
  
  
  F_23 <- 0.25 *(
    tr(sigma_matrix_inv %*% sigma2_2_deriv %*% 
         sigma_matrix_inv %*% a_deriv) - 
      B_sigma2_2 * B_a
    
  )
  
  F_24 <- 0.25 *(
    tr(sigma_matrix_inv %*% sigma2_2_deriv %*% 
         sigma_matrix_inv %*% rho_deriv) - 
      B_sigma2_2 * B_rho
    
  ) 
  
  #' # Fisher information third line
  #' 
  #' 
  #' 
  #' ## Fisher information for $\rho$
  
  F_33 <- 0.25 *(
    tr(sigma_matrix_inv %*% a_deriv %*% 
         sigma_matrix_inv %*% a_deriv) - 
      B_a^2
    
  )
  
  #' 
  #' ## Fisher information for $\rho$ and $a$
  #' 
  #' 
  #' 
  
  F_34 <- 0.25 *(
    tr(sigma_matrix_inv %*% a_deriv %*% 
         sigma_matrix_inv %*% rho_deriv) - 
      B_a * B_rho
    
  )
  
  
  #' # Fisher information fourth line
  #' 
  #' 
  
  F_44 <- 0.25 *(
    tr(sigma_matrix_inv %*% rho_deriv %*% 
         sigma_matrix_inv %*% rho_deriv) - 
      B_rho^2
    
  )
  
  # --------------- Assembling the expected values -----------
  
  fisher_info <- matrix(0,4,4)
  
  fisher_info[1,1] <- F_11
  fisher_info[1,2] <- F_12
  fisher_info[1,3] <- F_13
  fisher_info[1,4] <- F_14
  
  fisher_info[2,1] <- F_12
  fisher_info[2,2] <- F_22
  fisher_info[2,3] <- F_23
  fisher_info[2,4] <- F_24
  
  fisher_info[3,1] <- F_13
  fisher_info[3,2] <- F_23
  fisher_info[3,3] <- F_33
  fisher_info[3,4] <- F_34
  
  fisher_info[4,1] <- F_14
  fisher_info[4,2] <- F_24
  fisher_info[4,3] <- F_34
  fisher_info[4,4] <- F_44
  
  fisher_info
}

# ---- test -----

source("a_deriv_symmetry.R")
source("sigma2_1_deriv_symmetry.R")
source("sigma2_2_deriv_symmetry.R")
source("rho_deriv_symmetry.R")

set.seed(123)
coords_matrix <- matrix(runif(2*100), ncol = 2)
d <- dist(coords_matrix)
sigmas <- c(1,1)
nus <- c(0.5,0.5);nus[3] <- mean(nus)
as <- rep(2,3)
rho <- 0.5
theta_vec <- c(sigmas,as[3],rho)

M_dash_nu_1 <- BivMaternEstim:::matern_deriv(a = as[1], 
                                             nu = nus[1], 
                                             coords_matrix = coords_matrix)

M_dash_nu_2 <- BivMaternEstim:::matern_deriv(a = as[2], 
                                             nu = nus[2], 
                                             coords_matrix = coords_matrix)

M_dash_nu_3 <- BivMaternEstim:::matern_deriv(a = as[3], 
                                             nu = nus[3], 
                                             coords_matrix = coords_matrix)

create_sigma2_1_deriv(sigmas,as,rho,coords_matrix)
create_sigma2_2_deriv(sigmas,as,rho,coords_matrix)
create_a_deriv(M_dash_nu_1,M_dash_nu_2,M_dash_nu_3,sigmas,as,nus)
create_rho_deriv(sigmas,as,rho,coords_matrix)

debugonce(get_fisher_info)
get_fisher_info(theta_vec = theta_vec,nus = nus,coords_matrix = coords_matrix)
