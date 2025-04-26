sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("../../../old_code/")

load("../../For_15-05/gradient/saved_sim/biwmSim_and_idx.rda")

theta_true <- c(1,1.5,1,0.7)
nu_vec_true <- rep(1,2)

theta_start <- c(0.5,3,2,0.3)

res <- fit_biwm(obs_vec = obs_val,
         coords_matrix = coords,
         theta0 = theta_start,
         nus = nu_vec_true,
         verbosity = 6)
