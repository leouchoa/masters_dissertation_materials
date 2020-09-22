library(RandomFields)
library(BivMaternEstim)
library(ggplot2)
library(tidyr)
source("cov_with_new_loc.R")
source("matern_cov_wrapper_krig.R")


alr_inv <- function(gp_pred){
  
  aux <- exp(gp_pred) / (1 + sum(exp(gp_pred)))
  
  return(
    c(aux,1 - aux[1] - aux[2])
  )
  
}

set.seed(123) # doing again because of high code usage

n_points <- 40

coords_random <- matrix(runif(2 * (n_points) ,-1,1), ncol = 2)