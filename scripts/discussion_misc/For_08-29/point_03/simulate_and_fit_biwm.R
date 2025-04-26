simulate_and_fit_biwm <-
  function(
    true_param = c(1.5,0.7,1.3,0.4),
    initial_guess = c(0.5,1.2,4,-0.2),
    seed = 123,
    n_size = 40
  ){

  set.seed(seed)

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

fitbiwm_sim_summary <- function(simulation_obj){
  setNames(
    as.data.frame(
      rbind(
        true_param = round(resultado$true_param,5),
        initial_guess = round(simulation_obj$initial_guess,5),
        estimated_generic = round(simulation_obj$generic_guess$theta,5),
        estimated_informed = round(simulation_obj$informed_guess$theta,5)
      )
    ),
    c("sigma2_1","sigma2_2","a","rho")
  )
}

# ------- USAGE --------

#' Some points:
#'
#' 1. It has a LOT of bias for all parameters.
#' 3. I've changed 0.005 to 0.0005 and not much changed
#' 4. It has tendency to break as sample size increases
#'

resultado <- simulate_and_fit_biwm(n = 56,initial_guess = c(0.5,1.2,1.3,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)

resultado <- simulate_and_fit_biwm(n = 56,initial_guess = c(0.5,1.2,4.4,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)

resultado <- simulate_and_fit_biwm(n = 56,initial_guess = c(0.5,1.2,4.4,-0.2), seed = 123456, true_param = c(0.2,1.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)

resultado <- simulate_and_fit_biwm(n = 56,initial_guess = c(1.2,1.2,4.4,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)

resultado <- simulate_and_fit_biwm(n = 56,initial_guess = c(0.5,1.2,3,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)



resultado <- simulate_and_fit_biwm(n = 56,initial_guess = c(0.5,1.2,1.3,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)

resultado <- simulate_and_fit_biwm(n = 80,initial_guess = c(0.5,1.2,1.3,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)

resultado <- simulate_and_fit_biwm(n = 90,initial_guess = c(0.5,1.2,1.3,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)

resultado <- simulate_and_fit_biwm(n = 100,initial_guess = c(0.5,1.2,1.3,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)

resultado <- simulate_and_fit_biwm(n = 120,initial_guess = c(0.5,1.2,1.3,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)

resultado <- simulate_and_fit_biwm(n = 180,initial_guess = c(0.5,1.2,1.3,-0.2), seed = 123456, true_param = c(0.2,2.7,1.3,-0.4))
fitbiwm_sim_summary(resultado)
