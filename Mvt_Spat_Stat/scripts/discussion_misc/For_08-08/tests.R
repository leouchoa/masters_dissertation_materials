S_inv_Deriv <- cbind(rbind(B_11_sigma1,t(B_12_sigma1)),rbind(B_12_sigma1,B_22_sigma1))

all.equal(S_inv_Deriv,
          S_inv %*% Deriv_matrix_sigma2_1)

tr(S_inv %*% Deriv_matrix_sigma2_1)
tr(S_inv_Deriv)

plot(diag(S_inv %*% Deriv_matrix_sigma2_1))
points(diag(S_inv_Deriv),col = "red")

plot(diag(S_inv %*% Deriv_matrix_sigma2_1),
     diag(S_inv_Deriv))
abline(0,1)


S_inv_Deriv_rho <- cbind(rbind(B_11_rho,t(B_12_rho)),rbind(B_12_rho,B_22_rho))

all.equal(S_inv_Deriv_rho,
          S_inv %*% Deriv_matrix_rho)

S_inv_Deriv_a <- cbind(rbind(B_11_a,t(B_12_a)),rbind(B_12_a,B_22_a))

all.equal(S_inv_Deriv_a,
          S_inv %*% Deriv_matrix_a)

plot(diag(S_inv %*% Deriv_matrix_a),
     diag(S_inv_Deriv_a))
abline(0,1)
