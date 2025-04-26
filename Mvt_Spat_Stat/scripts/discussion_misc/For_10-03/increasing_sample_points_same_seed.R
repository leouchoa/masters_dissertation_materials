library(BivMaternEstim)
source("gen_initial_points.R")
source("simulate_and_estimate_biwm.R")

if(!dir.exists("increasing_sampled_points_same_seed")){
  dir.create("increasing_sampled_points_same_seed")
}

initial_point <- gen_initial_points(c(1,2,1.5,0.4))
the_seed <- sample.int(1e3, 1)


simulate_and_estimate_biwm(n_points = 40,
                           true_params = c(1,2,1.5,0.4),
                           initial_params = initial_point,
                           seed_number = the_seed
)

# ggsave("increasing_sampled_points_same_seed/biwm_sim_param_recover_sampled_size_40.pdf",width = 5,height = 8)

simulate_and_estimate_biwm(n_points = 80,
                           true_params = c(1,2,1.5,0.4),
                           initial_params = initial_point,
                           seed_number = the_seed
)

# ggsave("increasing_sampled_points_same_seed/biwm_sim_param_recover_sampled_size_80.pdf",width = 5,height = 8)

simulate_and_estimate_biwm(n_points = 100,
                           true_params = c(1,2,1.5,0.4),
                           initial_params = initial_point,
                           seed_number = the_seed
)


# ggsave("increasing_sampled_points_same_seed/biwm_sim_param_recover_sampled_size_100.pdf",width = 5,height = 8)

simulate_and_estimate_biwm(n_points = 110,
                           true_params = c(1,2,1.5,0.4),
                           initial_params = initial_point,
                           seed_number = the_seed
)

# ggsave("increasing_sampled_points_same_seed/biwm_sim_param_recover_sampled_size_110.pdf",width = 5,height = 8)
