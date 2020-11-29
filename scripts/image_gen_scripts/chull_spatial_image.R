library(ggplot2)
X <- matrix(stats::rnorm(2000), ncol = 2)
# chull(X)
# plot(X, cex = 0.5)
hpts <- chull(X)
hpts <- c(hpts, hpts[1])
# lines(X[hpts, ])

random_points_within_chull <- sample.int(nrow(X),4)
aux_df <- as.data.frame(X[random_points_within_chull, ])

chull_df <- as.data.frame(X[hpts,])

ggplot(chull_df,aes(V1,V2)) + 
  geom_polygon(
    alpha = 0.3
    ) + 
  theme_void() + 
  geom_point(
    data = aux_df[1:3,],
    aes(V1,V2), 
    size = 4
  ) +
  geom_point(
    data = aux_df[4,],
    aes(V1,V2),
    shape = "X", 
    size = 6
  )
