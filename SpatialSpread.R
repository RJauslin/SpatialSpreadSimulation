
pathIni <- getwd()
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

## ----setup, include=FALSE------------------------------------------------

#install.packages("Rcpp")
#install.packages("WaveSampling")
#install.packages("sampling")
#install.packages("BalancedSampling")
#install.packages("spatstat")
#install.packages("spsurvey")
#install.packages("sp")
#
#install.packages("knitr")
#install.packages("readr")
#install.packages("tikzDevice")
#
#install.packages("raster")
#install.packages("rgdal")
#install.packages("sf")
#install.packages("rgeos")
#install.packages("gridExtra")
#install.packages("ggvoronoi")
#install.packages("grid")
#install.packages("lattice")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("scico")
#
#install.packages("kableExtra")
#install.packages("magrittr")
#install.packages("plyr")

# LOAD PACKAGES
library(Rcpp)
library(WaveSampling)
library(sampling)
library(BalancedSampling)
library(spatstat)
library(spsurvey)
library(sp)

library(knitr)
library(readr)
library(tikzDevice)

library(raster)
library(rgdal)
library(sf)
library(rgeos)
library(gridExtra)
library(ggvoronoi)
library(grid)
library(lattice)
library(ggplot2)
library(ggrepel)
library(scico)

library(kableExtra)
library(magrittr)
library(plyr)


# theme plot of the paper
theme_wave <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family="sans",color = "black",size = 9),
      panel.spacing = unit(2, "lines"),
      # title
      plot.title = element_text(hjust = 0.5,size = 9),
      # axes
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      # legend
      legend.position="bottom",
      legend.title = element_text(size = 9,vjust = +1.0),
      legend.key.size = unit(0.3, "cm"),
      legend.key.width = unit(0.7,"cm") ,
      # background colors
      panel.background=element_blank(),
      panel.border=element_rect(colour = "black",fill = "transparent"),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      # keep edge black facet_wrap
      # strip.background = element_rect(fill="white"),
      strip.text =element_text(color = "black",size = 8)
      )
}



## ----GRTS,echo=FALSE-----------------------------------------------------
# grts function
GRTS = function(p,x,y){
  quiet <- function(x){
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }
  N = length(p)
  n = round(sum(p))
  index = 1:N
  s = rep(0,times=N)
  att = data.frame(x=x,y=y,mdcaty=p,ids=index);
  design = list(None=list(panel=c(Panel1=n), seltype="Continuous", caty.n=c("Caty 1"=n), over=0));
  res=quiet(spsurvey::grts(design, DesignID="Site", SiteBegin=1, type.frame="finite",
                     src.frame="att.frame", in.shape=NULL, sp.object=NULL, att.frame=att,
                     id=NULL, xcoord="x", ycoord="y", stratum=NULL, mdcaty="mdcaty", startlev=NULL,
                     maxlev=11, maxtry=1000, shift.grid=TRUE, do.sample=rep(TRUE, length(design)),
                     shapefile=FALSE, prjfilename=NULL, out.shape="sample"));

  s[res$ids]=1;
  s;
}


## ----distance,echo = FALSE,warning=FALSE,message = FALSE,results='hide',cache=TRUE----

# spatial configuration
X <- expand.grid(seq(1,3,0.05),seq(1,3,0.05))
Xp <- expand.grid(1:3,1:3)
N <- nrow(X)

# empty array for the distance
d <- array(rep(0,N*N),c(N,N))
d2 <-  array(rep(0,N*N),c(N,N))
d3 <-  array(rep(0,N*N),c(N,N))


# shift
perturb1 <- 1/12
perturb2 <- 1/4

# distance matrix
for(i in 1:N){
d[i,] <- distUnitk(as.matrix(X),k =i,tore = TRUE,toreBound = 3)
d2[i,] <- distUnitk(as.matrix(X),k =i,tore = FALSE,toreBound = 0)
}

# add shift
for(i in 1:nrow(X)){
tmp <- as.matrix(X)
tmp[i,1] <- tmp[i,1] + perturb1
tmp[i,2] <- tmp[i,2] + perturb2
d3[i,] <- distUnitk(tmp,k =i,tore = TRUE,toreBound = 3)
}

# data preparation for the plot
Xdat <- data.frame(x = X[,1],y = X[,2],v = d[1,],typ = rep("Tore",N))
Xdat <- rbind(Xdat,data.frame(x = X[,1],y = X[,2],v = d2[1,],typ = rep("Euclidean",N)))
Xdat <- rbind(Xdat,data.frame(x = X[,1],y = X[,2],v = d3[1,],typ = rep("Shifted Tore",N)))
Xdat$typ <- factor(Xdat$typ, levels = c("Euclidean", "Tore", "Shifted Tore"))
Xp <- data.frame(x = Xp[,1],y = Xp[,2])
perturb <- data.frame(x = perturb1+1,y = perturb2+1,v = 1,typ = "Shifted Tore")

# ggplot

# tikz(file = "distance.tex", width =  5.78851, height = 3,standAlone = FALSE)
p_distance <- ggplot(Xdat) +
geom_tile(aes(x = x,y = y,fill = v))+
geom_point(data = Xp,aes(x = x,y = y),size = 1,pch = 1,stroke = 0.5,colour = "grey50")+
geom_point(data = perturb,aes(x = x,y = y,fill = v),colour = "red",size = 1)+
facet_wrap(~typ)+
coord_fixed()+
scale_x_continuous(name="$x$",breaks = c(1,2,3) ,limits=c(0.95, 3.05)) +
scale_y_continuous(name="$y$",breaks = c(1,2,3) ,limits=c(0.95, 3.05)) +
labs(fill="$m^2$")+
scale_fill_scico(palette = "vik")+
  theme_wave()
print(p_distance)
# dev.off()


## ----stratification,echo = FALSE,warning=FALSE,message = FALSE,results='hide',cache = TRUE----

# spatial configuration
N <- 3*3
n <- 3
X <- as.matrix(expand.grid(seq(1,sqrt(N),1),seq(1,sqrt(N),1)))
pik <- rep(n/N,N)
D <- as(diag(1/pik),"sparseMatrix")

# W matrices and data preparation 
W1 <- round(wpik(X,pik,bound = 1,tore = TRUE,shift = FALSE,toreBound =sqrt(N)),2)
W2 <- round(wpik(X,pik,bound = 1,tore = FALSE,shift = FALSE,toreBound = 0),2)
W3 <- round(wpik(X,pik,bound = 1,tore = TRUE,shift = TRUE,toreBound = sqrt(N)),2)
dat1 <- summary(W1)
dat1$type <- rep("Tore",nrow(dat1))
dat2 <- summary(W2)
dat2$type <- rep("Euclidean",nrow(dat2))
dat3 <- summary(W3)
dat3$type <- rep("Shifted Tore",nrow(dat3))
dat <- rbind(dat1,dat2,dat3)
dat$type <- factor(dat$type,levels = c("Euclidean","Tore","Shifted Tore"))


# ggplot 

# tikz(file = "strat.tex", width = 5.78851, height = 3,standAlone = FALSE)
p_strat <- ggplot(data = dat) +
  geom_tile(aes(x =j,y = i, fill = as.factor(x)))+ scale_y_reverse()+
  scale_fill_grey(start = 0.8,end = 0.2,labels = c("1/6", "2/9", "1/3"))+
  coord_fixed()+
  facet_wrap(~type)+
  labs(fill = "$w_{kl}$")+
  theme_wave()+
  theme(axis.text.x=element_blank(),
       axis.text.y=element_blank(),
       axis.title.x=element_blank(),
       axis.title.y=element_blank())
print(p_strat)
# dev.off()


## ----Wsp2,echo = FALSE,warning=FALSE,results='hide',cache = TRUE---------

# spatial configuration

N_moran <- 250
n_moran <- 10
X_moran <- data.frame(x = runif(N_moran),y = runif(N_moran))
pik_moran <- rep(n_moran/N_moran,N_moran)

# W matrix
W_moran <- wpik(as.matrix(X_moran),pik_moran,tore = FALSE)

# convex hull for all strata into a data frame
st <- list()
for(i in 1:N_moran){
tmp <- X_moran[which(W_moran[i,]!=0),]
st[[i]] <- tmp[chull(tmp),]
st[[i]] <- cbind(st[[i]],n = i)
}
h <- data.frame()
for(i in 1:N_moran){
h <- rbind(h,st[[i]])
}

# sample selection 
s_moran <- wave(as.matrix(X_moran),pik_moran)

# convex hull of only selected units
index <- which(s_moran == 1)
dats <- data.frame()
for(i in 1:length(index)){
dats <- rbind(dats,h[which(h$n == index[i]),])
}


# plot strata
p_1 <- ggplot()+
geom_polygon(data = dats,aes(x = x,y = y,group = n),fill = "grey23",colour = "black",alpha = 0.4)+
geom_point(data = X_moran[which(s_moran == 1),],aes(x = x,y = y),pch = 16)+
geom_point(data = X_moran,aes(x = x,y = y),pch = 1,alpha = 0.4)+
geom_point(data = X_moran[s_moran == 1,],aes(x = x,y = y),pch = 16)+
coord_fixed()+
ggtitle("Intial strata")+
theme_wave()+
theme(axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank())


# matrix A
A_moran <- round(W_moran%*%as(diag(1/pik_moran),"sparseMatrix"),9)
dat_moran <- summary(A_moran)
dat_moran$type <- rep("Tore",nrow(dat_moran))

# plot matrix A
p_2 <-ggplot(data = dat_moran) +
geom_tile(aes(x =j,y = i, fill = as.factor(x)))+
scale_y_reverse()+
scale_fill_grey(start = 0.1,end = 0.1)+
coord_fixed()+
labs(fill = "$a_{kl}$")+
ggtitle("Stratification matrix")+
theme_wave()+
theme(axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position = "none")


# ggplot 

# tikz(file = "stratmat.tex", width = 5.78851, height = 3,standAlone = FALSE)
grid.arrange(p_1, p_2, nrow = 1)
# dev.off()


## ----voro1,echo=FALSE, warning=FALSE,message=FALSE,results='hide',cache = TRUE----

# spatial configuration and sample selection

vor <-  data.frame(x = runif(50),y = runif(50))
pik <- rep(20/50,50)
s1 <-  srswor(20,50)
s2 <- wave(as.matrix(vor),pik)

# data preparation
vor <- rbind(vor,vor)
name <- c(rep("Simple random sampling",50),
rep("Weakly associated vectors",50))
vor <- cbind(vor,name = name)
vor <- cbind(vor, s =  c(s1,s2))

# spatial balance
sb1 <- sb_vk(pik,as.matrix(vor[1:50,1:2]),s1)
sb2 <- sb_vk(pik,as.matrix(vor[1:50,1:2]),s2)
vor <- cbind(vor, v = c(sb1,sb2))

# grid for the voronoi ggplot
grid <- cbind(x = seq(-0.1,1.1,length.out = 100),y = rep(-0.1,100))
grid <- rbind(grid,cbind(x = rep(1.1,100), y = seq(-0.1,1.1,length.out = 100)))
grid <- rbind(grid,cbind(x =  seq(1.1,-0.1,length.out = 100), y = rep(1.1,100)))
grid <- rbind(grid,cbind(x =  rep(-0.1,100), y =seq(1.1,-0.1,length.out = 100)))
grid <- cbind(grid , group = rep(1,400))
grid <- as.data.frame(grid)

# tikz(file = "voro.tex", width = 5.78851, height = 3,standAlone = FALSE)
p_voro <- ggplot() +
  geom_voronoi(data = vor[vor$s == 1,],aes(x = x,y = y,fill = v),outline =grid,size = 0.1,alpha = 1,colour ="black")+
  geom_point(data = vor,aes(x = x,y = y),size = 1,shape = 1,alpha = 0.5)+
  geom_point(data = vor[vor$s == 1,] ,aes(x = x,y = y),size = 1,shape = 16,alpha = 1)+
  scale_fill_scico(palette = "vik",begin = 0.7,end = 0.0,limits = c(min(pik),max(sb1)))+
  facet_wrap(~name)+
  coord_fixed()+
  labs(fill = "$v_k$")+
  theme_wave()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.border=element_blank())
print(p_voro)
# dev.off()


## ----datageneration,echo=FALSE,cache=TRUE,warning=FALSE,message=FALSE----

# Data loading/generation--------

N <- 225
## CSR
# pp <- rpoispp(N)
# X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
load("./results/X_pp_16_10_2019.rda")
N_pp <- nrow(X_pp)

## Neyman-Scott
# nclust <- function(x0, y0, radius, n) {
# return(runifdisc(n, radius, centre=c(x0, y0)))
# }
# ns <- rNeymanScott(15, 1, nclust, radius=0.055, n=15)
# X_ns <- matrix(cbind(ns$x,ns$y),nrow = ns$n,ncol = 2)
load("./results/X_ns_16_10_2019.rda")
N_ns <- nrow(X_ns)

## grid
x <- seq(1,sqrt(N),1)
X_grid <- as.matrix(expand.grid(x,x))
N_grid <- N

# target sample size
n = 75

# pik
pik_pp <- rep(n/N_pp,N_pp)
pik_ns <- rep(n/N_ns,N_ns)
pik_grid <- rep(n/N,N)


# sample
s_pp <- wave(X_pp,pik_pp,tore = FALSE,shift = FALSE)
s_ns <- wave(X_ns,pik_ns,tore = FALSE,shift = FALSE)
s_grid <- wave(X_grid,pik_grid,tore = TRUE,shift =TRUE)



## ----artificialPlot,echo = FALSE, warning=FALSE,message=FALSE,results='hide',cache = TRUE----

# data preparation
dat1 <- data.frame(x = X_pp[,1],y = X_pp[,2],type = rep("Complete spatial randomness",nrow(X_pp)))
dat2 <- data.frame(x = X_ns[,1],y = X_ns[,2],type = rep("Neyman-Scott process",nrow(X_ns)))
dat3 <- data.frame(x = X_grid[,1],y = X_grid[,2],type = rep("Simple regular grid",nrow(X_grid)))
dat <- rbind(dat1,dat2,dat3)


# tikz(file = "artificial.tex", width = 5.78851, height = 2.1,standAlone = FALSE)
  p_simu <- ggplot()+
    geom_point(data = dat,
               aes(x = x,y = y,group = type),pch = 1,alpha = 1,size = 1,stroke = 0.5)+
    geom_point(data = dat1[which(s_pp == 1),],aes(x = x,y = y,group = type),pch = 16,size = 1)+
    geom_point(data = dat2[which(s_ns == 1),],aes(x = x,y = y,group = type),pch = 16,size = 1)+
    geom_point(data = dat3[which(s_grid == 1),],aes(x = x,y = y,group = type),pch = 16,size = 1)+
    facet_wrap(~type,scales = "free",ncol = 3) +
    theme_wave()+
    theme(axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      strip.text =element_text(color = "black",size = 7.5))
	print(p_simu)
# dev.off()


## ----SIMU_pp,echo=FALSE,message=FALSE,results='hide',cache = TRUE,warning=FALSE----

# simulation number and sample size
SIM <- 10000
n1 <- 25
n2 <- 50
n3 <- 100
n <- c(n1,n2,n3)

#LOAD data OR simualtion, the simulation is time-consuming so I saved the data in rds/csv files.

readRDS(file = "./results/pik_16_10_2019.rds")


## How the simulations are performed

# pik_pp <- list(rep(n1/N_pp,N_pp),
#                rep(n2/N_pp,N_pp),
#                rep(n3/N_pp,N_pp),
#                inclusionprobabilities(runif(N_pp),n1),
#                inclusionprobabilities(runif(N_pp),n2),
#                inclusionprobabilities(runif(N_pp),n3))
#
#
# pik_ns <- list(rep(n1/N_ns,N_ns),
#                rep(n2/N_ns,N_ns),
#                rep(n3/N_ns,N_ns),
#                inclusionprobabilities(runif(N_ns),n1),
#                inclusionprobabilities(runif(N_ns),n2),
#                inclusionprobabilities(runif(N_ns),n3))
#
# pik_grid <- list(rep(n1/N_grid,N_grid),
#                  rep(n2/N_grid,N_grid),
#                  rep(n3/N_grid,N_grid),
#                  inclusionprobabilities(runif(N_grid),n1),
#                  inclusionprobabilities(runif(N_grid),n2),
#                  inclusionprobabilities(runif(N_grid),n3))
#
#
#
#
#
# pik <- list(pik_pp,pik_ns,pik_grid)
# X <- list(X_pp,X_ns,X_grid)
# N <- list(N_pp,N_ns,N_grid)
#
#
#
# names(N) <- names(X) <- names(pik) <- c("pp","ns","grid")
#
# # prepare s and W
# s <-  rep(list(0),3)
# Ib1 <- Ib <- Sb <- W <- W1 <- rep(list(0),3)
# names(Ib1) <- names(Ib) <- names(Sb) <- names(W) <-  names(W1) <- names <-  c("pp","ns","grid")
# for(i in 1:3){
#   s[[i]] <-  W1[[i]]<- W[[i]]<- Ib[[i]] <- Ib1[[i]] <- Sb[[i]] <- rep(list(0),6)
#   names(s[[i]]) <- names(W[[i]]) <-names(W1[[i]]) <- names(Ib1[[i]]) <- names(Ib[[i]]) <- names(Sb[[i]]) <-c("n1","n2","n3","n1_un","n2_un","n3_un")
#
#   for(j in 1:length(s[[i]])){
#     Ib1[[i]][[j]] <- Ib[[i]][[j]] <- Sb[[i]][[j]] <- data.frame(matrix(rep(0,5*SIM),nrow = SIM,ncol =5))
#     colnames(Ib1[[i]][[j]]) <- colnames(Ib[[i]][[j]]) <- colnames(Sb[[i]][[j]]) <-  c("wave","lpm1","scps","srswor","grts")
#
#     s[[i]][[j]] <-  rep(list(0),5)
#     names(s[[i]][[j]]) <- c("wave","lpm1","scps","srswor","grts")
#     # pikt <- UPMEpiktildefrompik(pik[[i]][[j]])
#     # w_maxent <- pikt/(1-pikt)
#     # q[[i]][[j]] <- UPMEqfromw(w_maxent,round(sum(pik[[i]][[j]])))
#
#   }
# }
#
#
# for(tt in 1:SIM){
#
#   # calculate sample
#   for(i in 1:length(s)){
#     for(j in 1:length(s[[i]])){
#       #W
#       W[[i]][[j]] <- wpik(X[[i]],pik[[i]][[j]],tore = F,shift = F)
#       W[[i]][[j]]  <- drop0(W[[i]][[j]] - diag(diag(W[[i]][[j]])), tol = 0, is.Csparse = NA)
#
#       #W1
#       W1[[i]][[j]] <- wpikInv(X[[i]],pik[[i]][[j]],tore = F,shift = F,toreBound = 0)
#
#
#       #wave
#       s[[i]][[j]][[1]] <- as.vector(wave(X[[i]],pik[[i]][[j]],bound = 1,tore = F,shift = F))
#
#       #lpm1
#       s[[i]][[j]][[2]] <- rep(0,N[[i]])
#       s[[i]][[j]][[2]][lpm1(pik[[i]][[j]],X[[i]])] <- 1
#
#       #scps
#       s[[i]][[j]][[3]] <- rep(0,N[[i]])
#       s[[i]][[j]][[3]][scps(pik[[i]][[j]],X[[i]])] <- 1
#
#       #GRTS
#       s[[i]][[j]][[4]] <- GRTS(pik[[i]][[j]],X[[i]][,1],X[[i]][,2])
#
#       #srswor
#       if(j == 1 | j == 2 | j == 3){
#         s[[i]][[j]][[5]] <- srswor(sum(pik[[i]][[j]]),N[[i]])
#       }else{
#         s[[i]][[j]][[5]] <- UPmaxentropy(pik[[i]][[j]])
#         # s[[j]][[5]] <- UPMEsfromq(q[[j]])
#       }
#
#
#     }
#   }
#
#   # fill Ib and Sb
#   for(i in 1:length(s)){
#     for(j in 1:length(s[[i]])){
#       for(k in 1:length(s[[i]][[j]])){
#         Ib1[[i]][[j]][tt,k] <- IB(W1[[i]][[j]],s[[i]][[j]][[k]])
#         Ib[[i]][[j]][tt,k] <- IB(W[[i]][[j]],s[[i]][[j]][[k]])
#         Sb[[i]][[j]][tt,k] <- sb(pik[[i]][[j]],X[[i]],which(s[[i]][[j]][[k]] == 1))
#       }
#     }
#   }
#
# }
#
# save <- function(res,name){
#   for(i in 1:length(Ib)){
#     tmp <- do.call(rbind,res[[i]])
#     rownames(tmp) <- NULL
#     tmp$type <- rep(c("n1","n2","n3","n1_un","n2_un","n3_un"),each = SIM)
#     write.table(tmp,
#                 file = file.path("./results",paste0(name,names(res)[i],".csv")),
#                 row.names=FALSE, na="",col.names=TRUE, sep=",")
#   }
# }
# save(Ib,"Ib")
# save(Ib1,"Ib1")
# save(Sb,"Sb")



# LOAD data and preparation of the table.

Ib <- list()

Ibpp <- read_csv("./results/Ibpp_16_10_2019.csv")
Ibpp$type <- factor(Ibpp$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ibpp <- split(as.data.frame(Ibpp),f = Ibpp$type, drop = TRUE)
Ibpp <- lapply(Ibpp, function(x) { x["type"] <- NULL; x })
Ibns <- read_csv("./results/Ibns_16_10_2019.csv")
Ibns$type <- factor(Ibns$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ibns <- split(as.data.frame(Ibns),f = Ibns$type, drop = TRUE)
Ibns <- lapply(Ibns, function(x) { x["type"] <- NULL; x })
Ibgrid <- read_csv("./results/Ibgrid_16_10_2019.csv")
Ibgrid$type <- factor(Ibgrid$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ibgrid <- split(as.data.frame(Ibgrid),f = Ibgrid$type, drop = TRUE)
Ibgrid <- lapply(Ibgrid, function(x) { x["type"] <- NULL; x })

Ib[[1]] <- Ibpp
Ib[[2]] <- Ibns
Ib[[3]] <- Ibgrid

names(Ib) <- c("pp","ns","grid")


Ib1 <- list()

Ib1pp <- read_csv("./results/Ib1pp_16_10_2019.csv")
Ib1pp$type <- factor(Ib1pp$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib1pp <- split(as.data.frame(Ib1pp),f = Ib1pp$type, drop = TRUE)
Ib1pp <- lapply(Ib1pp, function(x) { x["type"] <- NULL; x })
Ib1ns <- read_csv("./results/Ib1ns_16_10_2019.csv")
Ib1ns$type <- factor(Ib1ns$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib1ns <- split(as.data.frame(Ib1ns),f = Ib1ns$type, drop = TRUE)
Ib1ns <- lapply(Ib1ns, function(x) { x["type"] <- NULL; x })
Ib1grid <- read_csv("./results/Ib1grid_16_10_2019.csv")
Ib1grid$type <- factor(Ib1grid$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib1grid <- split(as.data.frame(Ib1grid),f = Ib1grid$type, drop = TRUE)
Ib1grid <- lapply(Ib1grid, function(x) { x["type"] <- NULL; x })

Ib1[[1]] <- Ib1pp
Ib1[[2]] <- Ib1ns
Ib1[[3]] <- Ib1grid

names(Ib1) <- c("pp","ns","grid")

Sb <- list()

Sbpp <- read_csv("./results/Sbpp_16_10_2019.csv")
Sbpp$type <- factor(Sbpp$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Sbpp <- split(as.data.frame(Sbpp),f = Sbpp$type, drop = TRUE)
Sbpp <- lapply(Sbpp, function(x) { x["type"] <- NULL; x })
Sbns <- read_csv("./results/Sbns_16_10_2019.csv")
Sbns$type <- factor(Sbns$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Sbns <- split(as.data.frame(Sbns),f = Sbns$type, drop = TRUE)
Sbns <- lapply(Sbns, function(x) { x["type"] <- NULL; x })
Sbgrid <- read_csv("./results/Sbgrid_16_10_2019.csv")
Sbgrid$type <- factor(Sbgrid$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Sbgrid <- split(as.data.frame(Sbgrid),f = Sbgrid$type, drop = TRUE)
Sbgrid <- lapply(Sbgrid, function(x) { x["type"] <- NULL; x })

Sb[[1]] <- Sbpp
Sb[[2]] <- Sbns
Sb[[3]] <- Sbgrid

names(Sb) <- c("pp","ns","grid")


# means on the different data/design
for(i in 1:3){
  Ib1[[i]] <- lapply(Ib1[[i]],colMeans)
  Ib[[i]] <- lapply(Ib[[i]],colMeans)
  Sb[[i]] <- lapply(Sb[[i]],colMeans)
}
nom <- c(paste("$n =",n1,"$"),paste("$n =",n2,"$"),paste("$n =",n3,"$"))

for(i in 1:3){
  
  Ib1[[i]] <- cbind(data.frame(matrix(c(Ib1[[i]][[1]],Ib1[[i]][[2]],Ib1[[i]][[3]]),nrow = 3,byrow = T )),
      data.frame(matrix(c(Ib1[[i]][[4]],Ib1[[i]][[5]],Ib1[[i]][[6]]),nrow = 3,byrow = T )))
  Ib1[[i]] <- cbind(nom,Ib1[[i]])
  Ib1[[i]] <- cbind(rep(names(Ib1)[i],3),Ib1[[i]])

  
  Ib[[i]] <- cbind(data.frame(matrix(c(Ib[[i]][[1]],Ib[[i]][[2]],Ib[[i]][[3]]),nrow = 3,byrow = T )),
      data.frame(matrix(c(Ib[[i]][[4]],Ib[[i]][[5]],Ib[[i]][[6]]),nrow = 3,byrow = T )))
  Ib[[i]] <- cbind(nom,Ib[[i]])
  Ib[[i]] <- cbind(rep(names(Ib)[i],3),Ib[[i]])
 
  Sb[[i]] <- cbind(data.frame(matrix(c(Sb[[i]][[1]],Sb[[i]][[2]],Sb[[i]][[3]]),nrow = 3,byrow = T )),
      data.frame(matrix(c(Sb[[i]][[4]],Sb[[i]][[5]],Sb[[i]][[6]]),nrow = 3,byrow = T )))
  Sb[[i]] <- cbind(nom,Sb[[i]])
  Sb[[i]] <- cbind(rep(names(Sb)[i],3),Sb[[i]])
   
  colnames(Ib1[[i]]) <-colnames(Ib[[i]]) <- colnames(Sb[[i]]) <-
    c("design","n","wave","lpm1","scps","srswor","grts","wave","lpm1","scps","maxent","grts")

}




## ----tablepp,echo=FALSE--------------------------------------------------

datpp <- rbind(Ib1$pp,Ib$pp,Sb$pp)
datpp <- datpp[,-1]
colnames(datpp) <- c(" ","wave","lpm1","scps","grts","srswor","wave","lpm1","scps","grts","maxent")

kable(datpp, format = "latex",digits = 3, booktabs = T, caption = paste0("Spreading measures results based on ",SIM," simulations on the Complete spatial randomness dataset. The population size is equal to ",N_pp,"."),row.names = FALSE,escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 5, "Unequal probabilities" = 5),escape = F) %>%
 add_header_above(c(" " = 1,"Sampling design" = 10)) %>%
 group_rows("$I_{B_1}$",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("$I_B$",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("$B$",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9)


## ----table1,echo=FALSE---------------------------------------------------

datIb1 <- rbind(Ib1[[1]],Ib1[[2]],Ib1[[3]])[,2:ncol(Ib1[[1]])]
colnames(datIb1) <- c("$I_{B_1}$","wave","lpm1","scps","srswor","grts","wave","lpm1","scps","maxent","grts")

tableIb1 <-  kable(datIb1, format = "latex",digits = 3, booktabs = T, caption = paste0("Results of ",SIM," simualtions of Moran's $I_{B_1}$ spatial measure \\eqref{eq:wpik1}.","The three spatial configurations of the Section \\ref{sec:artificial} are taken into account."),escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 5, "Unequal probabilities" = 5),escape = F) %>%
 add_header_above(c(" " = 1,"Sampling design" = 10)) %>%
 group_rows("Complete spatial randomness",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Neyman-Scott process",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Simple regular grid",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9)


## ----table2,echo=FALSE---------------------------------------------------

datIb <- rbind(Ib[[1]],Ib[[2]],Ib[[3]])[,2:ncol(Ib[[1]])]
colnames(datIb) <- c("$I_{B}$","wave","lpm1","scps","srswor","grts","wave","lpm1","scps","maxent","grts")

tableIb <- kable(datIb, format = "latex",digits = 3, booktabs = T, caption = paste0("Results of ",SIM," simulations of Moran's $I_B$ spatial measure \\eqref{eq:wpik}.","The three spatial configurations of the Section \\ref{sec:artificial} are taken into account."),row.names = FALSE,escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 5, "Unequal probabilities" = 5)) %>%
 add_header_above(c(" " = 1,"Sampling design" = 10)) %>%
 group_rows("Complete spatial randomness",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Neyman-Scott process",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Simple regular grid",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9)


## ----table3,echo=FALSE---------------------------------------------------

datSb <- rbind(Sb[[1]],Sb[[2]],Sb[[3]])[,2:ncol(Sb[[1]])]
colnames(datSb) <- c("$B$","wave","lpm1","scps","srswor","grts","wave","lpm1","scps","maxent","grts")

tableSb <- kable(datSb, format = "latex",digits = 3, booktabs = T, caption = paste0("Results of ",SIM," simulations of the spatial balance measure $B$ \\eqref{eq:voro}.","The three spatial configurations of the Section \\ref{sec:artificial} are taken into account."),row.names = FALSE,escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 5, "Unequal probabilities" = 5)) %>%
 add_header_above(c(" " = 1,"Sampling design" = 10)) %>%
 group_rows("Complete spatial randomness",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Neyman-Scott process",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Simple regular grid",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9)


## ----meuse_plot,echo =FALSE----------------------------------------------

# data loading
data("meuse")
data("meuse.riv")
meuse.riv <- meuse.riv[which(meuse.riv[,2] < 334200 & meuse.riv[,2] > 329400),]

X <- scale(as.matrix(meuse[,1:2]))
pik <- inclusionprobabilities(meuse$copper,30)


## ----meusefig,echo = FALSE, warning=FALSE,message=FALSE,results='hide',cache = TRUE----

# transform data into a simple feature
meuse_sf <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")

# sample selection
s <- wave(X,pik)

# tikz(file = "meuse.tex", width =5, height = 4,standAlone = FALSE)
meuseplot <- ggplot()+
  geom_sf(data = meuse_sf,aes(size=copper),show.legend = 'point',shape = 1,stroke = 0.3)+
  geom_polygon(data = data.frame(x = meuse.riv[,1],y = meuse.riv[,2]),
               aes(x = x,y = y),
               fill = "lightskyblue2",
               colour= "grey50")+
  geom_point(data = meuse,
             aes(x = x,y = y,size = copper),
             shape = 1,
             stroke = 0.3)+
  geom_point(data = meuse[which(s == 1),],
             aes(x = x,y = y,size = copper),
             shape = 16)+
  labs(x = "Longitude",
       y = "Latitude",
       # title = "Meuse river",
       size = "Copper",
       caption = NULL)+
  scale_size(range = c(0.5, 3.5))+
   scale_y_continuous(breaks = c(50.96,50.97,50.98,50.99),
                      labels = c("$50.96^\\circ$N","$50.97^\\circ$N","$50.98^\\circ$N","$50.99^\\circ$N")) +
  scale_x_continuous(breaks = c(5.72,5.73,5.74,5.75,5.76),
                     labels = c("$5.72^\\circ$E","$5.73^\\circ$E","$5.74^\\circ$E","$5.75^\\circ$E","$5.76^\\circ$E"))+
  # coord_equal()+
  theme_wave()+
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey80"))
print(meuseplot)
# dev.off()


## ----spread_measure_meuse,echo = FALSE,message=FALSE,warning=FALSE,results='hide',cache = TRUE----
SIM <- 10000
n1 <- 15
n2 <- 30
n3 <- 75
n <- c(n1,n2,n3)
N_meuse <- nrow(X)
y <- meuse$cadmium*0.01


## LOAD data OR simualtion, the simulation is time-consuming so I saved the data in rds/csv files.
## How simulations are performed


# pik_meuse <- list(rep(n1/N_meuse,N_meuse),
#                rep(n2/N_meuse,N_meuse),
#                rep(n3/N_meuse,N_meuse),
#                inclusionprobabilities(meuse$copper,n1),
#                inclusionprobabilities(meuse$copper,n2),
#                inclusionprobabilities(meuse$copper,n3))
#
#
#
#
# s <- W1<- W <- Ib <- Ib1 <- Sb  <- y_hat <- sim <- V <- rep(list(0),6)
# names(s) <- names(W) <-names(W1) <- names(Ib1) <- names(Ib) <- names(Sb) <-names(y_hat)<- names(sim) <-names(V) <- c("n1","n2","n3","n1_un","n2_un","n3_un")
# for(j in 1:length(s)){
#   Ib1[[j]] <- Ib[[j]] <- Sb[[j]] <- sim[[j]] <- V[[j]] <- y_hat[[j]] <- data.frame(matrix(rep(0,5*SIM),nrow = SIM,ncol =5))
#   colnames(Ib1[[j]]) <- colnames(Ib[[j]]) <- colnames(Sb[[j]]) <- colnames(sim[[j]]) <- colnames(V[[j]]) <- colnames(y_hat[[j]]) <-  c("wave","lpm1","scps","grts","srswor")
#   s[[j]] <-  rep(list(0),5)
#   names(s[[j]]) <- c("wave","lpm1","scps","grts","srswor")
# }
#
#
#
# for(tt in 1:SIM){
#   for(j in 1:length(s)){
#     #W
#     W[[j]] <- wpik(X,pik_meuse[[j]],tore = F,shift = F)
#     W[[j]]  <- drop0(W[[j]] - diag(diag(W[[j]])), tol = 0, is.Csparse = NA)
#     #W1
#     W1[[j]] <- wpikInv(X,pik_meuse[[j]],tore = F,shift = F,toreBound = 0)
#     #wave
#     s[[j]][[1]] <- as.vector(wave(X,pik_meuse[[j]],bound = 1,tore = F,shift = F))
#
#     #lpm1
#     s[[j]][[2]] <- rep(0,N_meuse)
#     s[[j]][[2]][lpm1(pik_meuse[[j]],X)] <- 1
#
#     #scps
#     s[[j]][[3]] <- rep(0,N_meuse)
#     s[[j]][[3]][scps(pik_meuse[[j]],X)] <- 1
#
#     #grts
#     s[[j]][[4]] <- GRTS(pik_meuse[[j]],X[,1],X[,2])
#
#     #srswor
#     if(j == 1 | j == 2 | j == 3){
#       s[[j]][[5]] <- srswor(sum(pik_meuse[[j]]),N_meuse)
#     }else{
#       s[[j]][[5]] <- UPmaxentropy(pik_meuse[[j]])
#       # s[[j]][[5]] <- UPMEsfromq(q[[j]])
#     }
#   }
#
#   # fill Ib and Sb
#   for(j in 1:length(s)){
#     for(k in 1:length(s[[j]])){
#       Ib1[[j]][tt,k] <- IB(W1[[j]],s[[j]][[k]])
#       Ib[[j]][tt,k] <- IB(W[[j]],s[[j]][[k]])
#       Sb[[j]][tt,k] <- sb(pik_meuse[[j]],X,which(s[[j]][[k]] == 1))
#     }
#   }
#
#
#    for(j in 1:length(s)){
#     for(k in 1:length(s[[j]])){
#
#       y_hat[[j]][tt,k] <- sum(y[which(s[[j]][[k]]==1)]/pik_meuse[[j]][which(s[[j]][[k]] ==1)])
#       sim[[j]][tt,k] <- (sum(y[which(s[[j]][[k]]==1)]/pik_meuse[[j]][which(s[[j]][[k]] ==1)]) - sum(y))^2
#
#       if(names(s[[j]])[k] == "srswor"){
#         V[[j]][tt,k] <- varHAJ(y[which(s[[j]][[k]]==1)],pik_meuse[[j]][which(s[[j]][[k]]==1)],which(s[[j]][[k]]==1))
#       }else if(names(s[[j]])[k] == "wave"){
#         wtmp <- localmean.weight(x = X[which(s[[j]][[k]]==1),1],
#                                   y = X[which(s[[j]][[k]]==1),2],
#                                   prb = pik_meuse[[j]][which(s[[j]][[k]]==1)],
#                                   nbh = 2)
#         V[[j]][tt,k] <- localmean.var(y[which(s[[j]][[k]]==1)]/pik_meuse[[j]][which(s[[j]][[k]]==1)],wtmp)
#       }else{
#         wtmp <- localmean.weight(x = X[which(s[[j]][[k]]==1),1],
#                                   y = X[which(s[[j]][[k]]==1),2],
#                                   prb = pik_meuse[[j]][which(s[[j]][[k]]==1)],
#                                   nbh = 4)
#         V[[j]][tt,k] <- localmean.var(y[which(s[[j]][[k]]==1)]/pik_meuse[[j]][which(s[[j]][[k]]==1)],wtmp)
#       }
#     }
#   }
#
# }
#

# save <- function(res,name){
#     tmp <- do.call(rbind,res)
#     rownames(tmp) <- NULL
#     tmp$type <- rep(c("n1","n2","n3","n1_un","n2_un","n3_un"),each = SIM)
#     write.table(tmp,
#                 file = file.path("./results",paste0(name,".csv")),
#                 row.names=FALSE, na="",col.names=TRUE, sep=",")
# }
#
# save(Ib,"IbMeuse")
# save(Ib1,"Ib1Meuse")
# save(Sb,"SbMeuse")
# save(V,"varMeuse")
# save(sim,"VsimuMeuse")
# save(y_hat,"YhtMeuse")


Ib <- read_csv("./results/IbMeuse_16_10_2019.csv")
Ib$type <- factor(Ib$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib <- split(as.data.frame(Ib),f = Ib$type, drop = TRUE)
Ib <- lapply(Ib, function(x) { x["type"] <- NULL; x })

Ib1 <- read_csv("./results/Ib1Meuse_16_10_2019.csv")
Ib1$type <- factor(Ib1$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib1 <- split(as.data.frame(Ib1),f = Ib1$type, drop = TRUE)
Ib1 <- lapply(Ib1, function(x) { x["type"] <- NULL; x })

Sb <- read_csv("./results/SbMeuse_16_10_2019.csv")
Sb$type <- factor(Sb$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Sb <- split(as.data.frame(Sb),f = Sb$type, drop = TRUE)
Sb <- lapply(Sb, function(x) { x["type"] <- NULL; x })

V <- read_csv("./results/varMeuse_16_10_2019.csv")
V$type <- factor(V$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
V <- split(as.data.frame(V),f = V$type, drop = TRUE)
V <- lapply(V, function(x) { x["type"] <- NULL; x })

sim <- read_csv("./results/VsimuMeuse_16_10_2019.csv")
sim$type <- factor(sim$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
sim <- split(as.data.frame(sim),f = sim$type, drop = TRUE)
sim <- lapply(sim, function(x) { x["type"] <- NULL; x })

y_hat <- read_csv("./results/YhtMeuse_16_10_2019.csv")
y_hat$type <- factor(y_hat$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
y_hat <- split(as.data.frame(y_hat),f = y_hat$type, drop = TRUE)
y_hat <- lapply(y_hat, function(x) { x["type"] <- NULL; x })


Ib1 <- lapply(Ib1,colMeans)
Ib <- lapply(Ib,colMeans)
Sb <- lapply(Sb,colMeans)

nom <- c(paste("$n =",n1,"$"),paste("$n =",n2,"$"),paste("$n =",n3,"$"))

Ib1 <- cbind(data.frame(matrix(c(Ib1[[1]],Ib1[[2]],Ib1[[3]]),nrow = 3,byrow = T )),
      data.frame(matrix(c(Ib1[[4]],Ib1[[5]],Ib1[[6]]),nrow = 3,byrow = T )))
Ib1 <- cbind(nom,Ib1)

Ib<- cbind(data.frame(matrix(c(Ib[[1]],Ib[[2]],Ib[[3]]),nrow = 3,byrow = T )),
      data.frame(matrix(c(Ib[[4]],Ib[[5]],Ib[[6]]),nrow = 3,byrow = T )))
Ib <- cbind(nom,Ib)

Sb <- cbind(data.frame(matrix(c(Sb[[1]],Sb[[2]],Sb[[3]]),nrow = 3,byrow = T )),
      data.frame(matrix(c(Sb[[4]],Sb[[5]],Sb[[6]]),nrow = 3,byrow = T )))
Sb <- cbind(nom,Sb)


colnames(Ib1) <-colnames(Ib) <- colnames(Sb) <-  c(" ","wave","lpm1","scps","grts","srswor","wave","lpm1","scps","grts","maxent")


## ----spread_meuse_table,echo=FALSE,results='asis'------------------------

datmeuse <- rbind(Ib1,Ib,Sb)
# saveRDS(datvar, file = "C:/Users/jauslinr/switchdrive/SpreadSampling/draft/table/SPR.rds")

kable(datmeuse, format = "latex",digits = 3, booktabs = T, caption = paste0("Spreading measures results based on ",SIM," simulations on the Meuse dataset. The population size is equal to ",N_meuse,"."),row.names = FALSE,escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 5, "Unequal probabilities" = 5),escape = F) %>%
 add_header_above(c(" " = 1,"Sampling design" = 10)) %>%
 group_rows("$I_{B_1}$",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("$I_B$",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("$B$",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9)



## ----var_meuse,echo = FALSE,message=FALSE,warning=FALSE,results='hide'----

# Coverage rate
Up <- Low <- y_hat
for(l in 1:6){
  for(r in 1:5){
    Up[[l]][,r] <- y_hat[[l]][,r] + qnorm(0.5 + (95/100)/2)*sqrt(V[[l]][,r])
    Low[[l]][,r] <- y_hat[[l]][,r] - qnorm(0.5 + (95/100)/2)*sqrt(V[[l]][,r])
  }
}

coverage <- Up
for(l in 1:6){
  for(r in 1:5){
    coverage[[l]][,r] <- (Low[[l]][,r] < sum(y) &  sum(y) < Up[[l]][,r])
  }
}

coverage_m <- lapply(coverage,colMeans)
sim_m <- lapply(sim,colMeans)
V_m <- lapply(V,colMeans)

# Ration between averages of estimated and simu
aver_var <- sim_m
for(l in 1:6){
  for(r in 1:5){
    aver_var[[l]][r] <- V_m[[l]][r]/sim_m[[l]][r]
  }
}

sim <- lapply(sim,colMeans)
V <- lapply(V,colMeans)


nom <- c(paste("$n =",n1,"$"),paste("$n =",n2,"$"),paste("$n =",n3,"$"))




sim<- cbind(data.frame(matrix(c(sim[[1]],sim[[2]],sim[[3]]),nrow = 3,byrow = T )),
            data.frame(matrix(c(sim[[4]],sim[[5]],sim[[6]]),nrow = 3,byrow = T )))
sim <- cbind(nom,sim)

V <- cbind(data.frame(matrix(c(V[[1]],V[[2]],V[[3]]),nrow = 3,byrow = T )),
             data.frame(matrix(c(V[[4]],V[[5]],V[[6]]),nrow = 3,byrow = T )))
V <- cbind(nom,V)

coverage_m <- cbind(data.frame(matrix(c(coverage_m[[1]],coverage_m[[2]],coverage_m[[3]]),nrow = 3,byrow = T )),
              data.frame(matrix(c(coverage_m[[4]],coverage_m[[5]],coverage_m[[6]]),nrow = 3,byrow = T )))
coverage_m <- cbind(nom,coverage_m)

aver_var <- cbind(data.frame(matrix(c(aver_var[[1]],aver_var[[2]],aver_var[[3]]),nrow = 3,byrow = T )),
              data.frame(matrix(c(aver_var[[4]],aver_var[[5]],aver_var[[6]]),nrow = 3,byrow = T )))
aver_var <- cbind(nom,aver_var)

colnames(sim) <-colnames(V) <- colnames(coverage_m) <-colnames(aver_var) <-  c(" ","wave","lpm1","scps","grts","srswor","wave","lpm1","scps","grts","maxent")


## ----tablevar,echo=FALSE-------------------------------------------------

datvar <- rbind(sim,V,coverage_m,aver_var)
# saveRDS(datvar, file = "../draft/table/VAR.rds")

kable(datvar, format = "latex",digits = 3, booktabs = T, caption = paste0("Results of ",SIM," simulations on Meuse dataset. The population size is equal to ",N_meuse,"."," $v_{sim}$ is equal to the variance approximated by the simulations \\eqref{eq:varSIM}. $v$ depends on the sampling design.
For the srswor and maxent methods, we used the estimator $v_{HAJ}$ \\eqref{eq:varHAJ} while for the other sampling designs, we use $v_{LM}$ \\eqref{eq:varLM}. In the wave sampling design, the parameter for the neighbouring is set to three instead of four. Coverage rate of the 95\\% confidence intervals are computed as well as the ratio between averages of $v$ and $v_{SIM}$."),row.names = FALSE,escape = FALSE)%>%
  add_header_above(c(" " = 1,"Equal probabilities" = 5, "Unequal probabilities" = 5),escape = F) %>%
  add_header_above(c(" " = 1,"Sampling design" = 10)) %>%
  group_rows("$v_{SIM}$",1,3,escape = F,bold = F,latex_gap_space = "1em")%>%
  group_rows("$v$",4,6,escape = F,bold = F,latex_gap_space = "1em")%>%
  group_rows("Coverage of the 95% confidence interval",7,9,escape = T,bold = F,latex_gap_space = "1em") %>%
  group_rows("Ratio $v/v_{SIM}$",10,12,escape = F,bold = F,latex_gap_space = "1em") %>%
  kable_styling(font_size = 9)



## ----ptable1,echo=FALSE--------------------------------------------------
tableIb1

## ----ptable2,echo=FALSE--------------------------------------------------
tableIb

## ----ptable3,echo=FALSE--------------------------------------------------
tableSb


## ----vorocountexample,echo=FALSE,warning=FALSE,message=FALSE,results='hide',cache=TRUE----

# spatial  configuration
x <- seq(1,6,1)
X <- as.matrix(expand.grid(x,x))
pik <- rep(6/36,36)

# grid for the Voronoi polygons (ggvoronoi option)
grid <- cbind(x = seq(0.8,6.2,length.out = 100),y = rep(0.8,100))
grid <- rbind(grid,cbind(x = rep(6.2,100), y = seq(0.8,6.2,length.out = 100)))
grid <- rbind(grid,cbind(x =  seq(6.2,0.8,length.out = 100), y = rep(6.2,100)))
grid <- rbind(grid,cbind(x =  rep(0.8,100), y =seq(6.2,0.8,length.out = 100)))
grid <- cbind(grid , group = rep(1,400))
grid <- as.data.frame(grid)

s_cluster <- c(rep(1,6),rep(0,30))
s_spread <- c(rep(0,2),
              1,
              rep(0,2),
              1,
              rep(0,7),
              1,
              rep(0,2),
              1,
              rep(0,7),
              1,
              rep(0,2),
              1,
              rep(0,8))


W <- wpik(X,pik)
W <- drop0(W- diag(diag(W)))

v_spread <- sb_vk(pik,X,s_spread)
v_cluster <- sb_vk(pik,X,s_cluster)

dat <- as.data.frame(rbind(X,X))
colnames(dat) <- c("x","y")
dat$type <- rep(c("Spreaded","Clustered"),each = 36,times = 1)
dat$s <- c(s_spread,s_cluster)
dat$v <- c(v_spread,v_cluster)

# tikz(file = "vorocounter.tex", width = 5.78851, height = 3,standAlone = FALSE)
p_counter <- ggplot()+
  geom_voronoi(data = dat[dat$s == 1,],aes(x = x,y = y,fill = v),outline = grid ,size = 0.1,alpha = 1,colour ="black")+
  geom_point(data = dat,aes(x = x,y = y),shape = 1,size = 1,alpha = 0.5)+
  geom_point(data = dat[dat$s ==1,],aes(x = x,y = y),shape = 16,size = 1,alpha = 1)+
  scale_fill_scico(palette = "vik",begin = 0.7,end = 0.35,limits = c(min(pik),max(v_spread)))+
  labs(fill = "$v_k$")+
  facet_wrap(~type)+
  coord_fixed()+
  theme_wave()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.border=element_blank())
print(p_counter)
# dev.off()

sb_cluster <- round(sb(pik,X,which(s_cluster ==1)),3)
sb_spread <- round(sb(pik,X,which(s_spread ==1)),3)

IB_cluster <- round(IB(W,s_cluster),3)
IB_spread <- round(IB(W,s_spread),3)



## ----vorocountexample2,echo=FALSE,warning=FALSE,message=FALSE,results='hide'----

# spatial configuration 1
x <- seq(1,6,1)
X <- as.matrix(expand.grid(x,x))
pik <- rep(4/36,36)
s <- c(rep(0,14),
       1,1,
       rep(0,4),
       1,1,
       rep(0,14))


# spatial configuration 2
X2 <-data.frame(x = c(2.25,5.25,3.5,
                     3.5,3.4,3.6),
               y = c(2.25,2.25,5.25,
                     3.64,3.4,3.4))

for(i in 4:6){
  Y1 <- rnorm(2,0,0.04) + X2[i,1]
  Y2 <- rnorm(2,0,0.04) + X2[i,2]
  tmp <- data.frame(x = Y1, y = Y2)
  X2 <- rbind(X2,tmp)
}


for(i in 1:3){
    Y1 <- rnorm(8,0,0.2) + X2[i,1]
    Y2 <- rnorm(8,0,0.2) + X2[i,2]
    tmp <- data.frame(x = Y1, y = Y2)
    X2 <- rbind(X2,tmp)
}



X2 <- as.matrix(X2)
pik2 <- rep(3/36,36)
s2 <- c(rep(0,3),rep(1,3),rep(0,30))


# grid for the Voronoi polygons (ggvoronoi option)
grid <- cbind(x = seq(0.8,6.2,length.out = 100),y = rep(0.8,100))
grid <- rbind(grid,cbind(x = rep(6.2,100), y = seq(0.8,6.2,length.out = 100)))
grid <- rbind(grid,cbind(x =  seq(6.2,0.8,length.out = 100), y = rep(6.2,100)))
grid <- rbind(grid,cbind(x =  rep(0.8,100), y =seq(6.2,0.8,length.out = 100)))
grid <- cbind(grid , group = rep(1,400))
grid <- as.data.frame(grid)


# W matrices for the IB measure
W <- wpik(X,pik)
W <- drop0(W- diag(diag(W)))

W2 <- wpik(X2,pik2)
W2 <- drop0(W2 - diag(diag(W2)))

v2 <- sb_vk(pik2,X2,s2)
v <- sb_vk(pik,X,s)

dat <- as.data.frame(rbind(X,as.matrix(X2)))
colnames(dat) <- c("x","y")
dat$type <- rep(c("Example 1","Example 2"),each = 36,times = 1)
dat$s <- c(s,s2)
dat$v <- c(v,v2)


# ggplot 

# tikz(file = "vorocounter2.tex", width = 5.78851, height = 3,standAlone = FALSE)
p_counter <- ggplot()+
  geom_voronoi(data = dat[dat$s == 1,],aes(x = x,y = y,fill = v),outline = grid ,size = 0.1,alpha = 1,colour ="black")+
  geom_point(data = dat,aes(x = x,y = y),shape = 1,size = 1,alpha = 0.5)+
  geom_point(data = dat[dat$s ==1,],aes(x = x,y = y),shape = 16,size = 1,alpha = 1)+
  scale_fill_scico(palette = "vik",begin = 0.7,end = 0.35,limits = c(min(pik),1.5883333))+
  labs(fill = "$v_k$")+
  facet_wrap(~type)+
  coord_fixed()+
  theme_wave()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.border=element_blank())
print(p_counter)
# dev.off()

sb1 <- round(sb(pik,X,which(s ==1)),3)
sb2 <- round(sb(pik2,X2,which(s2 ==1)),3)

IB_1 <- round(IB(W,s),3)
IB_2 <- round(IB(W2,s2),3)


#

setwd(pathIni)


