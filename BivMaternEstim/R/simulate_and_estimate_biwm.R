#' Compositional Krigging Simulation for the Bivariate (Wittle-)Matern Model.
#'
#' @export
#'
#' @examples
#'
#' initial_point <- c(0.2756810, 0.7721986, 2.9531961, 0.9619193)
#' sim_result <- simulate_and_estimate_biwm(n_points = 80,
#' true_params = c(1,3,1,0.9),
#' initial_params = initial_point,
#' seed_number = sample.int(123123,1),
#' nus_vec = c(1,1.5)
#' )

simulate_and_estimate_biwm <- function(n_points,
                                       true_params,
                                       initial_params,
                                       nus_vec = c(0.5,0.5),
                                       seed_number = 123
                                       ){

  set.seed(seed_number)
  suppressPackageStartupMessages(library(RandomFields))
  suppressPackageStartupMessages(library(tidyr))

  true_sigmas <- true_params[1:2]
  true_a <- true_params[3]
  true_rho <- true_params[4]

  initial_sigmas <- initial_params[1:2]
  initial_a <- initial_params[3]
  initial_rho <- initial_params[4]

  alr_inv <- function(gp_pred_with_no_coords){


    alr_inv_low_level <- function(gp_pred){

      aux <- exp(gp_pred) / (1 + sum(exp(gp_pred)))

      return(
        c(aux,1 - aux[1] - aux[2])
      )

    }

    return(
      setNames(
        as.data.frame(
          t(
            apply(
              gp_pred_with_no_coords,
              1,
              alr_inv_low_level
            )
          )
        ),
        paste0("comp_0",1:3)
      )
    )

  }


  x <- seq(-1, 1, length.out = 40)
  model <- RMbiwm(nudiag = nus_vec, nured = 1, rhored = true_rho, cdiag = true_sigmas, s = rep(true_a,3))

  print(model)

  biwm_sim <- RFsimulate(model, x, x)
  my_grid <- expand.grid(x,x)
  biwm_sim_df <-
    data.frame(
      process_1 = biwm_sim$variable1,
      process_2 = biwm_sim$variable2,
      coord_x = my_grid$Var1,
      coord_y = my_grid$Var2
    )



  biwm_sim_df_comp <-
    cbind(
      alr_inv(biwm_sim_df[,1:2]),
      biwm_sim_df[,3:4]
      ) %>%
    setNames(
      c(paste0("comp_",1:3),"coord_x","coord_y")
    )

  true_compositions <- biwm_sim_df_comp %>%
    gather("key","Valor",1:3) %>%
    ggplot(aes(coord_x,coord_y, fill = Valor)) +
    geom_tile() +
    facet_wrap(~key, labeller =
                 labeller(key =c(comp_1 = "Composição 1",
                                 comp_2 = "Composição 2",
                                 comp_3 = "Composição 3"))) +
    labs(x = "", y = "") +
    scale_fill_viridis_c() +
    labs(
      title = "Processo Simulado"
    )


  sampled_biwm_sim_df <- biwm_sim_df[sample(nrow(biwm_sim_df),n_points),]

  generic_test <- fit_biwm(sampled_biwm_sim_df[,1:2],
                           sampled_biwm_sim_df[,3:4],
                           initial_params,
                           nus_vec
  )


  cat(
    "Estimated parameters are: \n \n",
    paste(
      round(
        generic_test$theta,
        3
      ),
      collapse =  "\n"
    ),
    "\n \n"
  )

  krig_grid <- as.matrix(expand.grid(x,x))


  sigma_hat <-
    sigma_assembler_biwm(
      sigmas = generic_test$theta[1:2],
      a = generic_test$theta[3]
      ,rho = generic_test$theta[4],
      coords_matrix = sampled_biwm_sim_df[,3:4],
      nus = nus_vec,
      combined = TRUE
    )

  sigma_hat_with_new_loc <-
    BivMaternEstim:::cov_with_new_loc(
      loc_new = krig_grid,
      loc_obs = as.matrix(sampled_biwm_sim_df[,3:4]),
      sigmas = generic_test$theta[1:2],
      a = generic_test$theta[3],
      rho = generic_test$theta[4],
      nus = nus_vec,
      combined = TRUE
    )

  z_prediction <-
    setNames(
      as.data.frame(
        cbind(
          t(
            apply(
              matrix(
                sigma_hat_with_new_loc %*%
                  solve(sigma_hat,
                        as.vector(
                          t(
                            apply(
                              sampled_biwm_sim_df[,1:2],
                              1,
                              function(x) x - generic_test$mu
                            )
                          )
                        )
                  ),
                ncol = 2
              ),
              1,
              function(x) x + generic_test$mu
            )
          ),
          krig_grid
        )
      ),
      c("process_01_pred","process_02_pred","coord_x","coord_y")
    )

  compositions_prediction <-
    setNames(
      cbind(
        alr_inv(z_prediction[,1:2]),
        krig_grid
      ),
      c("comp_01_pred","comp_02_pred","comp_03_pred" ,"coord_x","coord_y")
    )

  predicted_compositions <- compositions_prediction %>%
    gather("key","Valor",1:3) %>%
    ggplot(aes(coord_x,coord_y, fill = Valor)) +
    geom_tile() +
    facet_wrap(~key, labeller =
                 labeller(key =c(comp_01_pred = "Composição 1",
                                 comp_02_pred = "Composição 2",
                                 comp_03_pred = "Composição 3"))) +
    labs(x = "", y = "") +
    scale_fill_viridis_c() +
    labs(
      title = "Predição Composicional",
      subtitle = paste("Pontos Amostrados:", n_points)
    )


  gridExtra::grid.arrange(
    true_compositions,
    predicted_compositions
  )

  return(
    list(
      biwm_sim_df = biwm_sim_df,
      compositions_prediction = compositions_prediction
    )
  )
}




