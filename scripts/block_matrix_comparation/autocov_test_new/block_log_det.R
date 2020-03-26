block_log_det <- function(matrix_list){
  
  if(length(matrix_list) != 3){
    stop("sigma_assembler_biwm didn't return 3 matrices",call. = FALSE)
  }
  
  log(det(
    matrix_list[[1]] - 
      matrix_list[[3]] %*% 
      solve(matrix_list[[2]],matrix_list[[3]])
  )) +
    
    log(det(matrix_list[[2]]))
  
}

# --- Sanity Check ----

# source("sigma_assembler_biwm.R")
# source("matrices_assembler.R")
# source("matern_cov_wrapper.R")
# 
# distances <- dist(matrix(runif(30),ncol = 2))
# theta <- c(1.5,2,2,0.5)
# #
# test_matrix <- sigma_assembler_biwm(theta[1:2],theta[3],theta[4],c(0.5,0.5),distances)
# 
# log(det(matrices_assembler(M_1 = test_matrix$C_11,M_2 = test_matrix$C_22, M_12 = test_matrix$C_12)))
# 
# block_log_det(test_matrix)
# 
# # ok good
# 
# A = test_matrix$C_12
# B = test_matrix$C_22
# all.equal(A%*%B, B%*%A)
# all.equal(B,solve(A,B)%*%A)