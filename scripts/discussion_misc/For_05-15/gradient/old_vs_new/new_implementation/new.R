library(BivMaternEstim)

load("../saved_sim/biwmSim_and_idx.rda")

theta_true <- c(1,1.5,1,0.7)
nu_vec_true <- rep(1,2)

theta_start <- c(0.5,3,2,0.3)

res <- fit_biwm(obs_matrix = obs_val,
                coords_matrix = coords,
                theta0 = theta_start,
                nus = nu_vec_true,
                verbosity = 6)
