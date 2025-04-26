matern_deriv <- function(coords_matrix, a, nu){

  dist_vec <- as.vector(
    as.matrix(dist(coords_matrix))
  )
  # i'm following 'fields' package routine

  # Revise this to avoid using besselK explicitly, only use Matern
  dist_vec[dist_vec == 0] <- 1e-10 # is there a way to avoid this?

  first_term <- ( 2^(nu - 1) * dist_vec^(nu) * a^(nu - 1) ) /
    gamma(nu)

  snd_term <- 2 * nu * besselK(a * dist_vec,nu) -
    a * dist_vec * besselK(a * dist_vec, nu + 1)


  return(
    matrix(first_term*snd_term,
           ncol = nrow(coords_matrix),
           nrow = nrow(coords_matrix)
    )
  )

}
