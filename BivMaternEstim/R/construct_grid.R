#' Construct a grid
#'
#' @export
#'
#'
construct_grid <- function(sigma2_1_vec,
                           sigma2_2_vec,
                           a_vec,
                           rho_vec){

  setNames(
    expand.grid(sigma2_1_vec,
                sigma2_2_vec,
                a_vec,
                rho_vec),
    c("sigma2_1","sigma2_2","a","rho"))

}
