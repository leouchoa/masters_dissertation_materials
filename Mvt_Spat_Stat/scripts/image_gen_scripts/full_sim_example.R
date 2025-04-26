library(BivMaternEstim)
library(gstat)
library(gridExtra)


initial_point <- c(0.2756810, 0.7721986, 2.9531961, 0.9619193)
sim_result <-
  simulate_and_estimate_biwm(
    n_points = 80,
    true_params = c(1,3,1,0.9),
    initial_params = initial_point,
    seed_number = sample.int(123123,1),
    nus_vec = c(1,1.5),
    nug_vec = c(0,0) #c(0.5,0.5)
  )

save(sim_result,file = "scripts/image_gen_scripts/simulate_and_estimate_biwm_result.rda")

# load("scripts/image_gen_scripts/simulate_and_estimate_biwm_result.rda")

loc_new <- sim_result$compositions_prediction[,4:5]
coords <- sim_result$loc_used

pred_std_err_mat <-
  get_pred_std_err(
    loc_new = loc_new,
    loc_obs = coords,
    biwm_fit = sim_result$biwm_fit,
    nus = c(1,1.5),
    nug_vec = c(0,0)
    )

krig_vals <- BivMaternEstim:::alr(sim_result$compositions_prediction[,1:3])


sigma_bar_hat <- get_approx_alrInv_cond_covvar(
  krig_values = krig_vals,
  loc_new = loc_new,
  loc_obs = coords,
  biwm_fit = sim_result$biwm_fit,
  nus = c(1,1.5),
  nug_vec = c(0,0)
)


cmp_sim <- sim_result$biwm_sim_df[,1:2] %>%  
  BivMaternEstim:::alr_inv() %>% 
  cbind(sim_result$biwm_sim_df[,3:4]) %>% 
  gather("key","Valor",1:3) %>%
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller =
               labeller(key =c(comp_01 = "Composição 1",
                               comp_02 = "Composição 2",
                               comp_03 = "Composição 3"
                               )
               )
  ) +
  labs(x = "", y = "") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = "Processo Composicional Simulado"
  )


alr_comp_sim <- sim_result$biwm_sim_df %>% 
  gather("key","Valor",1:2) %>%
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller =
               labeller(key =c(process_1 = "Processo 1",
                               process_2 = "Processo 2"
                               )
                        )
             ) +
  labs(x = "", y = "") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = "Processo Composicional Transformado"
  )

cmp_sim_pred <- sim_result$compositions_prediction %>% 
  gather("key","Valor",1:3) %>%
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller =
               labeller(key =c(comp_01_pred = "Composição 1",
                               comp_02_pred = "Composição 2",
                               comp_03_pred = "Composição 3"))) +
  labs(x = "", y = "") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = "Predição Composicional"
  )

cmp_sim_pred_err_no_points <- sigma_bar_hat %>%
  gather("key","Valor",1:3) %>%
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller =
               labeller(key =c(comp_01_sd = "Composição 1 (E.P)",
                               comp_02_sd = "Composição 2 (E.P)",
                               comp_03_sd = "Composição 3 (E.P)"))) +
  labs(x = "", y = "") +
  scale_fill_viridis_c() +
  theme_minimal() +
  guides(color = FALSE) +
  labs(
    title = "Erro Padrão da Predição Composicional"
  )


cmp_sim_pred_err_with_points <- sigma_bar_hat %>%
  gather("key","Valor",1:3) %>%
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller =
               labeller(key =c(comp_01_sd = "Composição 1 (E.P)",
                               comp_02_sd = "Composição 2 (E.P)",
                               comp_03_sd = "Composição 3 (E.P)"))) +
  labs(x = "", y = "") +
  scale_fill_viridis_c() +
  theme_minimal() +
  geom_point(
    data = cbind(coords,Valor = 0.05),
    aes(coord_x,coord_y,col = "red")
    ) +
  guides(color = FALSE) +
  labs(
    title = "Erro Padrão da Predição Composicional"
  )


grid.arrange(cmp_sim_pred,cmp_sim_pred_err_with_points)

sampled_coords <- sim_result$loc_used

aux_biwm_sim_df <- 
  dplyr::inner_join(sampled_coords,sim_result$biwm_sim_df)

coordinates(aux_biwm_sim_df) <-  ~coord_x+coord_y

g <- gstat(NULL, 
           id = "process_1", 
           formula = process_1 ~ 1,
           data = aux_biwm_sim_df)
g <- gstat(g, 
           id = "process_2", 
           formula = process_2 ~ 1,
           data = aux_biwm_sim_df)

fct_labels <- c(process_1.process_2 = "Variograma Cruzado",
                process_1 = "Variograma Processo 1",
                process_2 = "Variograma Processo 2")

cross_vario <- variogram(g)
levels(cross_vario$id) <- 
  c("process_1","process_2","process_1.process_2")


ggplot(cross_vario[,c(2,3,6)], aes(dist,gamma)) + 
  geom_point() + 
  facet_wrap(~id, 
             labeller = labeller(id = fct_labels)
  ) + 
  theme_minimal() + 
  labs(
    x = "Distância",
    y = expression(2 * gamma[ij](h))
  )
