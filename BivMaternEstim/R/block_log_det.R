block_log_det <- function(matrix_list){

  if(length(matrix_list) != 3){
    stop("sigma_assembler_biwm didn't return 3 matrices",call. = FALSE)
  }

  log(det(
    matrix_list[[1]] -
      matrix_list[[3]] %*%
      solve(matrix_list[[2]],t(matrix_list[[3]]))
  )) +

    log(det(matrix_list[[2]]))

}

# --- Sanity Check ----
#
# library(fields)
# source("sigma_assembler_biwm.R")
# source("matern_cov_wrapper.R")
#
# set.seed(123)
# distances <- dist(matrix(runif(30),ncol = 2))
# theta <- c(1.5,2,2,0.5)
#
# test_matrix <- sigma_assembler_biwm(theta[1:2],theta[3],theta[4],c(0.5,0.5),distances)
#
# test_matrix_combined <- sigma_assembler_biwm(theta[1:2],theta[3],theta[4],c(0.5,0.5),distances,combined = TRUE)
#
# results_new_way <- block_log_det(test_matrix)
# results_old_way <- log(det(test_matrix_combined))
#
# round(results_new_way - results_old_way,10)
