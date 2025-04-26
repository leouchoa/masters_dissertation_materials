tr <- function(in_matrix){
  sum(diag(in_matrix))
}

quad_form <- function(y_left,in_matrix,y_right){
  t(y_left) %*% in_matrix %*% y_right
}