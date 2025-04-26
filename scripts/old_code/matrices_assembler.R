matrices_assembler <- function(M_1,M_2,M_12){
  cbind(rbind(M_1, t(M_12)), rbind(M_12, M_2))
}
