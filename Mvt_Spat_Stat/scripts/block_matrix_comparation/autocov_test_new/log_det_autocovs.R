source("sigma_assembler_biwm.R")
source("matern_cov_wrapper.R")
source("block_log_det.R")

get_log_dets <- function(theta,
         nus,
         dists){
  
  sigmas <- theta[1:2]
  as <- theta[3]
  rho <- theta[4]
  
  block_autocov_matrix <- sigma_assembler_biwm(sigmas,as,rho,nus,
                                      coords_matrix = dists)
  
  block_log_det(block_autocov_matrix)
}


# distances <- dist(matrix(runif(30),ncol = 2))

nus <- c(0.5,0.5)
sigma_1 <- seq(0.5,5,length.out = 20)
sigma_2 <- seq(0.5,5,length.out = 20)
as <- seq(0.5,2,length.out = 20)
rho = 0.5
malha <- expand.grid(sigma1 = sigma_1,sigma2 = sigma_2,a = as,rho = rho)

# ----- loop -----
  
do_it <- function(n_size_vec,save_folder,n_seeds=3){
  
  for (j in seq_along(n_size_vec)) {
    
    det_out_list <- vector("list",n_seeds)
    
    for(k in 1:n_seeds){
      
      det_out <- vector("double",nrow(malha))
      set.seed(k)
      distances <- dist(matrix(runif(2*n_size_vec[j]),ncol = 2))
      
      for(i in seq_len(nrow(malha))){
        theta_in <- c(sigmas = c(malha[i,1], malha[i,2]), 
                      a = malha[i,3], 
                      rho = malha[i,4])
        
        det_out[i] <- get_log_dets(theta = theta_in,nus = nus,dists = distances)
      }
      
      det_out_list[[k]] <- det_out
    }
    
    save(det_out_list, 
         file = paste0(save_folder,"n_size_",n_size_vec[j],".rda"))
    
  }
}

#------------------- end

#---- Usage ----

# n_size_vec <- seq(110,300,by=10)
# save_folder <- "autocovs_dets/"
# 
# do_it(n_size_vec,save_folder,2)



arqs <- dir("autocovs_dets",full.names = TRUE)

sum_isInf_matrix <- matrix(0,length(arqs),3)

for(i in seq_along(arqs)){
  
  load(arqs[i])
  
  sum_isInf_matrix[i,1] <- as.integer(
    stringr::str_extract(arqs[i],"\\d+"))
  
  sum_isInf_matrix[i,2:ncol(sum_isInf_matrix)] <- colSums(apply(
    do.call(cbind,det_out_list),2,is.infinite))
  
}

sum_isInf_matrix <- as.data.frame(sum_isInf_matrix)
sum_isInf_matrix <- sum_isInf_matrix[order(sum_isInf_matrix$V1),]

