krig_cv_fn <- function(one_out_idx){
  
  fit_dts <- soil_dts@data[-one_out_idx,]
  
  fit_dts_coords <- coordinates(soil_dts)[-one_out_idx,]
  
  krig_new_point <- coordinates(soil_dts)[one_out_idx, , drop = FALSE]
  
  soil_dts_fit <- fit_biwm(
    obs_matrix = fit_dts,
    coords_matrix = fit_dts_coords,
    theta0 = initial_point,
    nus = nus_vec,
    nug_vec = nug_vec
  )  
  
  #otherwise, cant invert matrix with just 1 data point......
  krig_new_point <- coordinates(soil_dts)[c(one_out_idx,1), , drop = FALSE]
  
  loocv_pred <- compositional_biwm_krig(
    biwm_fit = soil_dts_fit,
    krig_locations = krig_new_point,
    fit_locations = fit_dts_coords,
    obs_matrix = fit_dts,
    nus = nus_vec,
    nug_vec
  )
  
  loocv_pred[1,]
}

cv_to_labels <- function(prop_dts){
  
  
  class_label <- as.data.frame(TT.points.in.classes(prop_dts*100,class.sys= "SiBCS13.TT"))
  
  
  
  apply(class_label,1,function(x){
    aux_term <- colnames(class_label)[which(as.logical(x))];
    ifelse(length(aux_term) != 1,paste0(aux_term,collapse = "/"),aux_term)
  })
  
  
}


library(BivMaternEstim)
library(gstat)
library(sp)
library(ggplot2)
library(tidyr)
library(ggmap)
library(viridis)
library(soiltexture)

KRIG_GRID_LEN <- 110

alr_transformed_with_locations_UNIQUE <- read.csv("data/alr_transformed_with_locations_UNIQUE.csv", row.names=NULL)

soil_dts <- alr_transformed_with_locations_UNIQUE[,1:4]

coordinates(soil_dts) <-  ~coord_x + coord_y


initial_point <- c(1,1,0.5,0)
nus_vec <- c(1,1)
nug_vec <- c(0.5,0.5)

loocv_krig_pred <- as.data.frame(matrix(0,ncol = 5,nrow = nrow(soil_dts)))

for (i in 1:nrow(loocv_krig_pred)) {
  
  loocv_krig_pred[i,] <- krig_cv_fn(i)
}

loocv_krig_pred <- setNames(loocv_krig_pred,
  c("SAND","CLAY","SILT","coord_x","coord_y")
)

save(loocv_krig_pred,file = "loocv_krig_pred.rda")



krig_labels <- cv_to_labels(loocv_krig_pred)
dts_labels <- cv_to_labels(
  setNames(alr_transformed_with_locations_UNIQUE[,8:10],
           c("SAND","CLAY","SILT"))
  )

table(krig_labels,dts_labels)
(12+45)/nrow(alr_transformed_with_locations_UNIQUE)