library(spdep)

rnd_coords <- matrix(runif(2*10),ncol = 2)

star_point <- rbind(c(0.25,0.25),c(0.5,0.5))

plot(rnd_coords)
points(star_point,pch = 4)

knn_res <- knearneigh(rnd_coords,3)

plot(knn2nb(knn_res),rnd_coords)

nn_mat <- nb2mat(knn2nb(knn_res))

rnd_coords_new <- rbind(matrix(runif(2*10),ncol = 2),star_point)

knn_res_new <- knearneigh(rnd_coords_new,3)

plot(knn2nb(knn_res_new),rnd_coords_new)

nn_mat_new <- nb2mat(knn2nb(knn_res_new))

idxs <- seq(nrow(rnd_coords_new) - nrow(star_point) +1 ,nrow(rnd_coords_new))

aux <- apply(nn_mat_new[idxs,] != 0,1,which)

for(i in seq_len(ncol(aux))){
  print(rnd_coords[aux[,i],])
}


# gotta use loops in the beginning. otherwise if the newest points are close together, we cant use them in the prediction


