biwm_plot <- function(distances,
                      sigma_1,
                      sigma_2,
                      a,
                      rho,
                      nus = c(0.5,0.5,0.5),
                      titulo=NULL,
                      ...){
  
  
  
  if(length(nus) != 3){
    stop("length(nus) != 3")
  }
  
  gen_form <- function(nu_i){
    
    fields::Matern(d = distances,
                         alpha = a,
                         nu = nu_i) 
  }
  
  dts <- data.frame(
    distances,
    C_11 = sigma_1*gen_form(nus[1]),
    C_12 = rho*sqrt(sigma_1)*sqrt(sigma_2)*gen_form(nus[3]),
    C_21 = rho*sqrt(sigma_1)*sqrt(sigma_2)*gen_form(nus[3]),
    C_22 = sigma_2*gen_form(nus[2])) %>%
    bind_rows(mutate(.,distances = -distances))
  
  library(magrittr)  # well yeah i'm lazy
  
  dts %>% 
    gather("key","val",-distances) %>% 
    ggplot(aes(distances,val)) +
    geom_line() +
    facet_wrap(~key, scales = "free_y") +
    labs(y = "", title = titulo, x= "Distancia")
  
}
