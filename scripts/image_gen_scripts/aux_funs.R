#' # Aux functions

#' ## biwm_plot
#' plots the supposed bivariate covariance function for many


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
  
  dts %>% 
    gather("key","val",-distances) %>% 
    ggplot(aes(distances,val)) +
    geom_line() +
    facet_wrap(~key, scales = "free_y") +
    labs(y = "", title = titulo, x= "Distancia")
  
}



#' ## is_valid_cov
#' test whether is a valid covariance function or not, given a value for rho and a vector for smoothing parameters.
#' nu_vec must be bivariate

is_valid_cov <- function(rho,nu_vec = c(0.5, 1)){
  test_res <- sqrt(prod(nu_vec))/mean(nu_vec)
  
  cat(
    paste("|rho| =",abs(rho),"is less than",test_res),"\n \n",
    "So is",abs(rho) <= sqrt(prod(nu_vec))/mean(nu_vec),
    "that this is a valid cov.fun"
    
  )
}


#' ## biwm_sim_rho
#' Simulate biwm processes from RandomFields for different values of \rho and some fixed params. It uses `map_df` from `purrr` package with a vector, here named rho_vec.
#' 
#' example: res <- map_df(rho_vec,biwm_sim_rho)
#' 
#' \textbf{important:} the `print_model = TRUE` to print the model used after each simulation is to be passed latter on but it is work in progress.

biwm_sim_rho <- function(rho_param,print_model = TRUE){
  model <- RMbiwm(nudiag=c(0.5, 0.5), nured=1, rhored=rho_param, cdiag=c(1, 1), s=c(1, 1, 1))
  
  sim <- RFsimulate(model, x, y)
  
  if(print_model){
    print(model)  
  }
  
  invisible(setNames(
    data.frame(sim,coordinates(sim)),
    c("Processo 01","Processo 02","coord_x","coord_y")
  ))
  
}

#' ## biwm_sim_many_rho_plot
#' 
#' Plots simulated biwm processes for various values of \rho and other parameter values given by `biwm_sim_rho`.

biwm_sim_many_rho_plot <- function(rho_vec=c(-0.8,-0.2,0.2,0.8),...){
  
  res <- map_df(rho_vec,biwm_sim_rho)
  
  # next loop creates the rho_val vector used to identify the simulation values to the associated rhos
  
  for(i in seq_along(rho_vec)){
    ifelse(i == 1,
           rho_val <- rep(rho_vec[i],nrow(res)/length(rho_vec)),
           rho_val <- append(rho_val,rep(rho_vec[i],nrow(res)/length(rho_vec))))
  }
  
  res$rho <- rho_val
  
  res %>% 
    gather("key","Valor",1:2) %>% 
    ggplot(aes(coord_x,coord_y, fill = Valor)) +
    geom_tile() +
    facet_wrap(key~rho,nrow=2) + 
    labs(x = "", y = "",
         title = "Simulacao Matern Bivariada para Varios valores de rho") +
    scale_fill_viridis_c() + coord_fixed()
  
}

#' ## biwm_sim_a
#' Simulate biwm processes from RandomFields for different values of \rho and some fixed params. It uses `map_df` from `purrr` package with a vector, here named rho_vec.
#' 
#' example: res <- map_df(rho_vec,biwm_sim_a)
#' 
#' \textbf{important:} the `print_model = TRUE` to print the model used after each simulation is to be passed latter on but it is work in progress.

biwm_sim_a <- function(a,print_model = TRUE){
  
  s_arg <- rep(1/a,3)
  model <- RMbiwm(nudiag=c(0.5, 0.5), nured=1, rhored=0.2, cdiag=c(1,1), s=s_arg)
  
  sim <- RFsimulate(model, x, y)
  
  if(print_model){
    print(model)
  }
  
  invisible(setNames(
    data.frame(sim,coordinates(sim)),
    c("Processo 01","Processo 02","coord_x","coord_y")
  ))
  
}


#' ## biwm_sim_many_a_plot
#' 
#' Plots simulated biwm processes for various values of \rho and other parameter values given by `biwm_sim_a`.

biwm_sim_many_a_plot <- function(a_vec = round(1 / seq(1,sqrt(800),length.out = 4), 4)){
  
  res <- map_df(a_vec,biwm_sim_a)
  
  # next loop creates the a_val vector used to identify the simulation values to the associated as
  
  for(i in seq_along(a_vec)){
    ifelse(i == 1,
           a_val <- rep(a_vec[i],nrow(res)/length(a_vec)),
           a_val <- append(a_val,rep(a_vec[i],nrow(res)/length(a_vec))))
  }
  
  res$a <- paste0("a=",a_val,"(s=",round(1/a_val,3),")")
  
  res %>% 
    gather("key","Valor",1:2) %>% 
    ggplot(aes(coord_x,coord_y, fill = Valor)) +
    geom_tile() +
    facet_wrap(key~a,nrow=2) + 
    labs(x = "", y = "",
         title = "Simulacao Matern Bivariada para Varios valores de a", subtitle = paste0(round(sqrt(800),4),"é a diagonal da caixa")) +
    scale_fill_viridis_c() + coord_fixed()
  
}


#' ## biwm_sim_nu
#' Simulate biwm processes from RandomFields for different values of $\nu_2$ and some fixed params. It uses `map_df` from `purrr` package with a vector, here named nu_2_vec.
#' 
#' example: res <- map_df(nu_2_vec,biwm_sim_nu)
#' 
#' \textbf{important:} the `print_model = TRUE` to print the model used after each simulation is to be passed latter on but it is work in progress.
#' 
#' 
biwm_sim_nu <- function(nu_2_vec,print_model = TRUE){
  
  model <- RMbiwm(nudiag=c(0.5, nu_2_vec), nured=1, rhored=0.2, cdiag=c(1, 1), s=c(1, 1, 1))
  
  sim <- RFsimulate(model, x, y)
  
  if(print_model){
    print(model)  
  }
  
  invisible(setNames(
    data.frame(sim,coordinates(sim)),
    c("Processo 01","Processo 02","coord_x","coord_y")
  ))
  
}
#' 
#' 
#' ## biwm_sim_many_nu_plot
#' 
#' Plots simulated biwm processes for various values of \rho and other parameter values given by `biwm_sim_nu`.
#' 
#' 
biwm_sim_many_nu_plot <- function(nu_2_vec,...){
  
  res <- map_df(nu_2_vec,biwm_sim_nu)
  
  # next loop creates the nu_2_val vector used to identify the simulation values to the associated rhos
  
  for(i in seq_along(nu_2_vec)){
    ifelse(i == 1,
           nu_2_val <- rep(nu_2_vec[i],nrow(res)/length(nu_2_vec)),
           nu_2_val <- append(nu_2_val,rep(nu_2_vec[i],nrow(res)/length(nu_2_vec))))
  }
  
  res$nu_2 <- paste("nu_2 = ", nu_2_val)
  res$nu <- paste("nu_1 = ", 0.5)
  
  res %>% 
    gather('key',"val",1:2) %>% 
    ggplot(aes(coord_x,coord_y,fill=val)) +
    geom_tile() +
    facet_wrap(key ~ nu_2 + nu) +
    labs(x = "", y = "",
         title = "Simulacao Matern Bivariada para Varios valores de nu") +
    scale_fill_viridis_c() + coord_fixed()
  
}
