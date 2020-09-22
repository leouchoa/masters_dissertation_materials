library(BivMaternEstim)

alr_inv <- function(gp_pred){
  
  aux <- exp(gp_pred) / (1 + sum(exp(gp_pred)))
  
  return(
    c(aux,1 - aux[1] - aux[2])
  )
  
}

set.seed(3)
n <- 40
coords <- matrix(runif(2*n), ncol = 2)
temp <- rnorm(2*n)
S <- sigma_assembler_biwm(sigmas = c(1, 1), a = 2, rho = 0.5,
                          nus = c(0.5, 0.5), coords_matrix = coords, combined = TRUE)
log_cd <- matrix(temp%*%chol(S) + rep(c(1,2), each = n), ncol = 2)

generic_test <- fit_biwm(log_cd, coords, c(.5, .5, 4, .6), c(0.5, 0.5))


x_grid <- seq(0,1,length.out = 10)
krig_grid <- as.matrix(expand.grid(x_grid,x_grid))
# krig_grid <- setNames(expand.grid(x_grid,x_grid),NULL)

# points(krig_grid[,1],krig_grid[,2],pch = "s")


sigma_hat <- 
  sigma_assembler_biwm(
    sigmas = generic_test$theta[1:2],
    a = generic_test$theta[3]
    ,rho = generic_test$theta[4],
    coords_matrix = coords,
    nus = c(0.5,0.5),
    combined = TRUE
  )

# workaround while it's not in the package
# source("../../../BivMaternEstim/R/matern_cov_wrapper.R")

sigma_hat_with_new_loc <- 
  cov_with_new_loc(
    loc_new = krig_grid,
    loc_obs = coords,
    sigmas = generic_test$theta[1:2],
    a = generic_test$theta[3],
    rho = generic_test$theta[4],
    nus = c(0.5,0.5),
    combined = TRUE
  )

z_prediction <- 
  setNames(
    as.data.frame(
      cbind(
        matrix(
          generic_test$mu + 
            sigma_hat_with_new_loc %*% 
            solve(sigma_hat,
                  as.vector(log_cd) - as.vector(t(apply(log_cd, 1, function(x) x - generic_test$mu)))
            ),
          ncol = 2
        ),
        krig_grid
      )
    ),
    c("process_01_pred","process_02_pred","coord_x","coord_y")
  )

compositions_prediction <- 
  setNames(
    cbind(
      as.data.frame(
        t(apply(as.data.frame(z_prediction[,1:2]),
                1,
                alr_inv)
        )
      ),
      krig_grid
    ),
    c("comp_01_pred","comp_02_pred","comp_03_pred" ,"coord_x","coord_y")
)

compositions_prediction %>% 
  gather("key","Valor",1:3) %>% 
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller = 
               labeller(key =c(comp_01_pred = "Composição 1",
                               comp_02_pred = "Composição 2",
                               comp_03_pred = "Composição 3"))) + 
  labs(x = "", y = "") +
  scale_fill_viridis_c()

  
  