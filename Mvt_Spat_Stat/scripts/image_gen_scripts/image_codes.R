#funcoes para crição de gráficos etc

# library(viridis)
library(tidyverse)
library(patchwork)
library(ggforce)
library(RandomFields)
library(gstat)
library(sp)
# library(ggtern) #please don't load this
theme_set(theme_minimal())
options(ggplot2.continuous.colour="viridis")
source("aux_funs.R")

# ----fig 1-----

a = data.frame(x = 1:2,y = 3:4)

ggplot(a,aes(x,y)) + 
  geom_point(size=3) + 
  scale_x_continuous(limits = c(-2,4)) + 
  scale_y_continuous(limits = c(-1,6)) + 
  geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="both", type = "closed")) + theme_void() + labs(x="",y="")

ggsave("Documents/dissertacao/texto_dissertacao/Diss_Leonardo/imagens/2pts_stationarity_v2.pdf")


# ----fig 2 -----

ggplot(a,aes(x,y)) + 
  geom_point(size=3) + 
  scale_x_continuous(limits = c(0,4)) + 
  scale_y_continuous(limits = c(1,6)) + 
  geom_line(color = "red") + 
  theme_minimal() + labs(x="",y="") + 
  theme(
    axis.text.x=element_blank(),
    axis.text.y = element_blank()
    )

ggsave("Documents/dissertacao/texto_dissertacao/Diss_Leonardo/imagens/2pts_isotropy.pdf")

# ----Drawing circles -----


circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# dat <- circleFun(c(1,3),1,npoints = 100)

ggplot(a,aes(x,y)) +
  geom_point(size=3) + 
  coord_fixed()+
  geom_line(color = "red") + 
  theme_minimal() + labs(x="",y="") + 
  theme(axis.text.x=element_blank(),axis.text.y = element_blank()) +
  geom_path(data = circleFun(c(1,3),1,npoints = 100),aes(x,y)) +
  geom_path(data = circleFun(c(1,3),1.5,npoints = 100),aes(x,y)) +
  geom_path(data = circleFun(c(1,3),0.5,npoints = 100),aes(x,y)) +
  geom_path(data = circleFun(c(1,3),2,npoints = 100),aes(x,y)) 

ggsave("imagens/2pts_isotropy.pdf")

ggplot(a,aes(x,y)) +
  geom_point(size=3) + 
  coord_fixed()+
  geom_line(color = "red") + 
  theme_minimal() + 
  labs(x="",y="") + 
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank()
        ) +
  geom_ellipse(aes(x0 = 1, y0 = 3, a = 1, b = 0.5, angle = pi / 4)) +
  geom_ellipse(aes(x0 = 1, y0 = 3, a = 0.8, b = 0.3, angle = pi / 4)) + 
  geom_ellipse(aes(x0 = 1, y0 = 3, a = 0.6, b = 0.15, angle = pi / 4))
ggsave("imagens/2pts_anisotropy.pdf")


p1 <- ggplot(a,aes(x,y)) +
  geom_point(size=3) + 
  coord_fixed()+
  geom_line(color = "red") + 
  theme_minimal() + labs(x="",y="") + 
  theme(axis.text.x=element_blank(),axis.text.y = element_blank()) +
  geom_path(data = circleFun(c(1,3),1,npoints = 100),aes(x,y)) +
  geom_path(data = circleFun(c(1,3),1.5,npoints = 100),aes(x,y)) +
  geom_path(data = circleFun(c(1,3),0.5,npoints = 100),aes(x,y)) +
  geom_path(data = circleFun(c(1,3),2,npoints = 100),aes(x,y)) +
  ggtitle("Isotropia")

p2 <- ggplot(a,aes(x,y)) +
  geom_point(size=3) + 
  coord_fixed()+
  geom_line(color = "red") + 
  theme_minimal() + labs(x="",y="") + 
  theme(axis.text.x=element_blank(),axis.text.y = element_blank()) +
  geom_ellipse(aes(x0 = 1, y0 = 3, a = 1, b = 0.5, angle = pi / 4)) +
  geom_ellipse(aes(x0 = 1, y0 = 3, a = 0.8, b = 0.3, angle = pi / 4)) + 
  geom_ellipse(aes(x0 = 1, y0 = 3, a = 0.6, b = 0.15, angle = pi / 4)) + 
  geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="both", type = "closed")) +
  ggtitle("Anisotropia")

p1 + p2

ggsave("imagens/2pts_isotropy&_anisotropy.pdf")


# ---- RF simulation -----

set.seed(123)

x <- seq(0, 1, len=100)
model <- RMwhittle(nu=1, Aniso=matrix(nc=2, c(1.5, 3, -3, 4)))
plot(model, dim=2, xlim=c(-1,1))
z <- RFsimulate(model=model, x, x)
plot(z)


# ---- Brownian Motion -------


t <- seq(0,2,length.out = 30)  # time
sigma2 <- 0.01
## first, simulate a set of random deviates
x <- rnorm(n = length(t) - 1, sd = sqrt(sig2))
## now compute their cumulative sum
x <- c(0, cumsum(x))
bmotion <- data.frame(time = t, bm = x)
ggplot(bmotion,aes(time,bm)) + 
  geom_line() +
  labs(x = "tempo", y = "B(tempo)")

ggsave("imagens/brownian_motion.pdf")

# ----- Gaussian Process ---

set.seed(123) # doing again because of high code usage

dists <- seq(0,5,0.1)
biwm_plot(dists,1,3,1,0.8,c(0.5,1.5,1.0))

example("simulate_and_estimate_biwm")

biwm_sim_df <- sim_result$biwm_sim_df

# ----- Gaussian Process Subset Estimation and Variogram ------

cbind(
  BivMaternEstim:::alr_inv(biwm_sim_df[,1:2]),
  biwm_sim_df[,3:4]
) %>% 
  gather("key","Valor",1:3) %>% 
  ggplot(aes(coord_x,coord_y, fill = Valor)) +
  geom_tile() +
  facet_wrap(~key, labeller = 
               labeller(key =c(comp_01 = "Composição 1",
                               comp_02 = "Composição 2",
                               comp_03 = "Composição 3"))) + 
  labs(x = "", y = "") +
  scale_fill_viridis_c()

idx <- sample(nrow(biwm_sim_df),80)

save(biwm_sim_df,idx,file = "biwm_sim_&_idx.rda")

biwm_sim_df_sampled <- biwm_sim_df[idx,]

biwm_sim_df_sampled %>% 
  gather("key","val",-(coord_x:coord_y)) %>% 
  ggplot(aes(coord_x,coord_y,color = val)) + 
  geom_point() +
  facet_wrap(~key)

# --- Variogram ---

coordinates(biwm_sim_df_sampled) <-  ~coord_x+coord_y

g <- gstat(NULL, 
           id = "process_1", 
           formula = process_1 ~ 1,
           data = biwm_sim_df_sampled)
g <- gstat(g, 
           id = "process_2", 
           formula = process_2 ~ 1,
           data = biwm_sim_df_sampled)

fct_labels <- c(process_1.process_2 = "Variograma Cruzado",
                process_1 = "Variograma Processo 1",
                process_2 = "Variograma Processo 2")

cross_vario <- variogram(g)
ggplot(cross_vario[,c(2,3,6)], aes(dist,gamma)) + 
  geom_point() + 
  facet_wrap(~id, 
             labeller = labeller(id = fct_labels)
             )
    # theme_minimal()


# ----- Apendix Simulations -----

set.seed(123)


x <- y <- seq(-10, 10, 0.1)

biwm_sim_many_rho_plot()
biwm_sim_many_a_plot()
biwm_sim_many_nus_plot() #not viable for the moment
biwm_sim_many_sigmas_plot() #not viable for the moment