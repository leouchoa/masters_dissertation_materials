# ---- Load -----

load('sim_res_kinda_same_seed.rda')
true_param <- c(1.5,0.7,1.3,0.4)

# ---- Generic Guess -----

par(mfrow = c(2,2), pty = "s")

hist(res_generic[1,],xlab = "");abline(v = true_param[1], col = "red",lwd = 5)
hist(res_generic[2,],xlab = "");abline(v = true_param[2], col = "red",lwd = 5)
hist(log(res_generic[3,]),xlab = "");abline(v = log(true_param[3]), col = "red",lwd = 5)
hist(res_generic[4,],xlab = "");abline(v = true_param[4], col = "red",lwd = 5)

# ---- Informed Guess -----

par(mfrow = c(2,2), pty = "s")

hist(res_informed[1,],xlab = "");abline(v = true_param[1], col = "red",lwd = 5)
hist(res_informed[2,],xlab = "");abline(v = true_param[2], col = "red",lwd = 5)
hist(res_informed[3,],xlab = "");abline(v = true_param[3], col = "red",lwd = 5)
hist(res_informed[4,],xlab = "");abline(v = true_param[4], col = "red",lwd = 5)