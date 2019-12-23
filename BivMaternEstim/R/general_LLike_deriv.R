general_LLike_deriv <- function(cov_matrix,
                                grad_matrix,
                                obs_vec){

  y <- solve(cov_matrix, obs_vec) # Need to center on mu

  -0.5 * (
    sum(diag(
      solve(cov_matrix, grad_matrix)
    )) -
      as.double(t(y) %*% grad_matrix %*% y)
  )

}
