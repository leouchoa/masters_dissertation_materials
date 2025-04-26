library(BivMaternEstim)

library(ggplot2)

true_theta <- c(1,1,2,0.5)
initial_pts <- BivMaternEstim:::gen_initial_points(true_theta)

theta_hat_hist_res_plot <- theta_hat_hist(
  true_theta,
  initial_pts,
  n_points = 80,
  nus_vec = c(0.5, 0.5),
  nug_vec = c(0, 0),
  n_replicates = 500,
  seed_number = sample.int(12312312, 1)
)

ggsave(filename = "../scripts/image_gen_scripts/theta_hat_hist_res_plot.pdf",
       width = 8.27,
       height = 11.69)

save(theta_hat_hist_res_plot,file = "../scripts/image_gen_scripts/theta_hat_hist_res_plot.rda")
