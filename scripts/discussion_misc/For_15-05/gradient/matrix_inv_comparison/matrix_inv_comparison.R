#' Guilherme aqui está um script pra comparar a inversão das matrizes. 
#' 
#' A única jogada aqui é que nas funções `sigma_assembler_biwm` e `block_inv` existe um argumento chamado `combined`. Ele serve para, ao invés de retornar uma lista com matrizes-blocos invertidas, retornar a inversa na forma de data.frame. Ou seja, se você olhar o apêndice B, `combined = FALSE` devolve cada elemento de $\Sigma^{-1}(\textbf{h})$ em (B.2). Por outro lado `combined = TRUE` devolve a matriz $\Sigma^{-1}(\textbf{h})$ de (B.2) como um data.frame.
#' 
#' #' 
#' 
#' 
#' 

suppressPackageStartupMessages(library(BivMaternEstim))

distances <- dist(matrix(runif(30),ncol = 2))
theta <- c(1.5,2,2,0.5)
nu_vec <- c(0.5,0.5)

test_matrix <- sigma_assembler_biwm(sigmas = theta[1:2],
                                    a = theta[3],
                                    rho = theta[4],
                                    nus = nu_vec,
                                    coords_matrix = distances)

test_matrix_combined <- sigma_assembler_biwm(sigmas = theta[1:2],
                                    a = theta[3],
                                    rho = theta[4],
                                    nus = nu_vec,
                                    coords_matrix = distances,
                                    combined = TRUE)



block_inv_test_matrix <- BivMaternEstim:::block_inv(test_matrix,
                                              combined = TRUE)

inv_test_matrix_combined <- solve(test_matrix_combined)

par(pty = "s",mfrow = c(1,2))
image(block_inv_test_matrix,main = "block inverse")
image(inv_test_matrix_combined,main = "R's solve")