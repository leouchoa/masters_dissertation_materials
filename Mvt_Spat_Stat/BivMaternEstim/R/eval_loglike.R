#' Eval it
#'
#' @export
#'
#'
eval_loglike <- function(theta_vec,log_cd, coords){

  LLike_biwm(theta = theta_vec,
             nus = c(0.5,0.5),
             mu = colMeans(log_cd),
             coords_matrix = coords,
             obs_matrix = log_cd)

}
