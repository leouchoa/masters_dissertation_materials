sigma_normal_llike <- function(sigma_sqrt){

  res <- rep(0,length(sigma_sqrt))

  for(i in seq_along(sigma_sqrt)){
    res[i] <- sum(log(dnorm(
      seq(-2,2,length.out = 50),
      sd=sigma_sqrt[i]
      )))
  }

  res
}

new_grid <- seq(0.1,5,length.out = 50)

sigma_sqrt_llike <- sigma_normal_llike(new_grid)

plot(new_grid,sigma_sqrt_llike)
