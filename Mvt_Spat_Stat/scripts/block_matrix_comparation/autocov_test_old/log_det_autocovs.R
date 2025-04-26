sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("source_files/")

get_dets <- function(theta,
         nus,
         coords_matrix){
  
  sigmas <- theta[1:2]
  as <- theta[3]
  rho <- theta[4]
  
  autocov_matrix <- sigma_assembler_biwm(sigmas = sigmas,
                                         a = as,
                                         rho = rho,
                                         nus = nus,
                                         coords_matrix = coords_matrix)
  
  det(autocov_matrix)
}


sigma_1 <- seq(0.5,5,length.out = 20)
sigma_2 <- seq(0.5,5,length.out = 20)
as <- seq(0.5,2,length.out = 20)
rho = 0.5
malha <- expand.grid(sigma1 = sigma_1,sigma2 = sigma_2,a = as,rho = rho)



do_it <- function(n_size_vec,save_folder,n_seeds=3){
  
  for (j in seq_along(n_size_vec)) {
    
    det_out_list <- vector("list",n_seeds)
    
    for(k in 1:n_seeds){
      
      det_out <- vector("double",nrow(malha))
      set.seed(k)
      coords <- matrix(runif(2*n_size_vec[j]),ncol = 2)
      
      for(i in 1:nrow(malha)){
        
        theta_in <- c(
          sigmas = c(malha[i,1], malha[i,2]), 
          a = malha[i,3], 
          rho = malha[i,4])
        
        det_out[i] <- get_dets(theta_in,c(0.5,0.5),coords)
        
      }
      
      det_out_list[[k]] <- det_out
    }
    
    save(det_out_list, 
         file = paste0(save_folder,"n_size_",n_size_vec[j],".rda"))
    
  }
}

#---- Usage ----

# n_size_vec <- seq(110,300,by=10)
# save_folder <- "autocovs_dets/"
# 
# do_it(n_size_vec,save_folder,3)

arqs <- dir("autocovs_dets",full.names = TRUE)

sum_log_matrix <- matrix(0,length(arqs),4)

sum_isInf_matrix <- matrix(0,length(arqs),4)

for(i in seq_along(arqs)){
  load(arqs[i])
  
  sum_isInf_matrix[i,1] <- as.integer(
    stringr::str_extract(arqs[i],"\\d+"))
  
  sum_isInf_matrix[i,2:4] <- colSums(apply(
    log(do.call(cbind,det_out_list)),2,is.infinite))
}

sum_isInf_matrix <- as.data.frame(sum_isInf_matrix)
sum_isInf_matrix <- sum_isInf_matrix[order(sum_isInf_matrix$V1),]


cbind(sum_isInf_matrix[,1],
      apply(sum_isInf_matrix[,2:4],2,function(x){x / 8000 * 100}))
