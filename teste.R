library(BivMaternEstim)

set.seed(1)
n <- 75
coords <- matrix(runif(2*n), ncol = 2)
temp <- rnorm(2*n)
S <- sigma_assembler_biwm(sigmas = c(1, 1), a = 2, rho = 0.5, nus = c(0.5, 0.5), coords_matrix = coords)
log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)

informed_test <- fit_biwm(log_cd, coords, c(1, 1, 2, .5), c(0.5, 0.5))
informed_test
generic_test <- fit_biwm(log_cd, coords, c(.5, .5, 4, .6), c(0.5, 0.5))
generic_test

B <- 100
sampleT <- matrix(NA_real_, ncol = 4, nrow = B)
sampleM <- matrix(NA_real_, ncol = 2, nrow = B)
for(i in 1:B){
  print(i)
  temp <- rnorm(2*n)
  log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)
  informed_test <- try({
    fit_biwm(log_cd, coords, c(1, 1, 2, .5), c(0.5, 0.5))
  })
  if(class(informed_test) == "try-error") next()
  sampleM[i,] <- informed_test$mu
  sampleT[i,] <- informed_test$theta
}

sampleM <- data.frame(mux = sampleM[,1], muy = sampleM[,2])
sampleT <- data.frame(sigma1 = sampleT[,1], sigma2 = sampleT[,2], a = sampleT[,3], rho = sampleT[,4])

require(ggplot2)
require(dplyr)
require(tidyr)
gather(sampleM) %>% ggplot(aes(x = key, y = value)) + geom_boxplot()
gather(sampleT) %>% ggplot(aes(x = key, y = value)) + geom_boxplot() + facet_wrap(~key, scales = "free_y")




optim(theta0,
      fn = LLike_biwm,
      gr = LLike_biwm_reduced_grad,
      method = "L-BFGS-B",
      lower = c(0.001, 0.001, 0.001, -min(sqrt(nus[1]*nus[2])/mean(nus), 0.999)),
      upper = c(Inf, Inf, Inf, min(sqrt(nus[1]*nus[2])/mean(nus), 0.999)),
      control = list(fnscale = -1, trace = 6),
      nus = nus,  mu = mu, # We always use Z = vec(\mathbf{Y} - boldsymbol\mu), maybe center it outside function?
      coords_matrix = coords_matrix,
      obs_vec = obs_vec)
# Olhar mvtnorm
