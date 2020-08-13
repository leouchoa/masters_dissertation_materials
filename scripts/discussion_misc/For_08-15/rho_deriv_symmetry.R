set.seed(123)
coords <- matrix(runif(2*100), ncol = 2)
d <- dist(coords)
sigmas <- c(1,1)
nus <- c(0.5,0.5);nus[3] <- mean(nus)
a <- 2
rho <- 0.5

D_11_rho <- matrix(0,nrow(coords),nrow(coords))

D_21_rho <- sqrt(sigmas[2]/ sigmas[1]) * 
  BivMaternEstim:::matern_cov_wrapper(d,
                                      a = a,
                                      nu = nus[3]
                                      )

D_22_rho <- matrix(0,nrow(coords),nrow(coords))

Deriv_matrix_rho <- cbind(rbind(D_11_rho, D_21_rho), rbind(t(D_21_rho), D_22_rho))

range(Deriv_matrix_rho - t(Deriv_matrix_rho))
