block_inv <- function(matrix_list,combined = FALSE){

  V_1 <- solve(matrix_list[[2]],matrix_list[[3]])

  V_2 <- matrix_list[[1]] - matrix_list[[3]] %*% V_1

  C_11_st <- solve(V_2)

  C_12_st <- - tcrossprod(C_11_st,V_1)

  # C_21_st <- t(C_12_star)

  C_22_st <- solve(matrix_list[[2]]) - V_1 %*% C_12_st

  out <- list(C_11_star = C_11_st,
              C_12_star = C_12_st,
              # C_21_star = t(C_12_st),
              C_22_star = C_22_st)

  if(combined){
    return(
      cbind(rbind(C_11_st, t(C_12_st)), rbind(C_12_st, C_22_st))
    )
  }else return(out)
}


# --- Sanity Check ----

# sourceDir <- function(path, trace = TRUE, ...) {
#   for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
#     if(trace) cat(nm,":")
#     source(file.path(path, nm), ...)
#     if(trace) cat("\n")
#   }
# }
#
# sourceDir(".")

# source("sigma_assembler_biwm.R")
# source("matrices_assembler.R")
# source("matern_cov_wrapper.R")
#
# distances <- dist(matrix(runif(30),ncol = 2))
# theta <- c(1.5,2,2,0.5)
#
# test_matrix <- sigma_assembler_biwm(theta[1:2],theta[3],theta[4],c(0.5,0.5),distances)
#
#
# block_inv(test_matrix)
#
# aux_assembler_inv <- function(matrix_list){
#   cbind(
#     rbind(matrix_list[[1]], t(matrix_list[[2]])),
#           rbind(matrix_list[[2]], matrix_list[[3]])
#     )
# }
#
# aux_assembler_autocov <- function(matrix_list_2){
#   cbind(
#     rbind(matrix_list_2[[1]], matrix_list_2[[3]]),
#     rbind(matrix_list_2[[3]], matrix_list_2[[2]])
#   )
# }
#
# a <- aux_assembler_inv(block_inv(test_matrix))
# b <- solve(aux_assembler_autocov(test_matrix))
# par(pty = "s",mfrow = c(1,2))
# image(t(a),main = "a")
# image(t(b),main = "b")
#
# all.equal(a,b)
#
#
# colSums(a == b)
# apply(a == b,2,which) # to-do latter
