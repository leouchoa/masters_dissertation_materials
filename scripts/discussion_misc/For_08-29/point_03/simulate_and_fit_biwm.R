simulate_and_fit_biwm <- function(true_param = c(1.5,0.7,1.3,0.4),
                        initial_guess = c(0.5,1.2,4,-0.2),
                        seed = 123,
                        n_size = 40){

  set.seed(seed)
  n_size <- 50

  coords <- matrix(runif(2*n_size), ncol = 2)

  temp <- rnorm(2*n_size)

  S <- sigma_assembler_biwm(sigmas = true_param[1:2], a = true_param[3], rho = true_param[4], nus = c(0.5, 0.5), coords_matrix = coords, combined = TRUE)

  log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n_size), ncol = 2)

  informed_test <- fit_biwm(log_cd, coords, true_param, c(0.5, 0.5))


  generic_test <- fit_biwm(log_cd, coords, initial_guess, c(0.5, 0.5))

  return(
    list(
      generic_guess = generic_test,
      informed_guess = informed_test,
      true_param = true_param,
      initial_guess = initial_guess
      )
    )
}

# ------- USAGE --------

#' Some points:
#'
#' 1. It has a LOT of bias for both $\sigma^2_1, \sigma^2_1$ and $a$.
#' 2. For $\rho$ it looks like the bias less than the others.
#' 3. I've changed 0.005 to 0.0005 and not much changed

resultado <- simulate_and_fit_biwm(n = 56,initial_guess = c(0.5,1.2,4.4,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
resultado$generic_guess
