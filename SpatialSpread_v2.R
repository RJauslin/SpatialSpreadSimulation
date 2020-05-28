
# These lines are here to set the working directory
# equal to the full path where SpatialSpread.R is saved.
#
# It allows that relative paths './results/' are working.
# You can remove these lines or specify the path where './results/' apprears.

pathIni <- getwd()
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)


## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------

# install.packages("Rcpp")
# # install.packages("WaveSampling")
# install.packages("sampling")
# install.packages("BalancedSampling")
# install.packages("spatstat")
# install.packages("spsurvey")
# install.packages("sp")
# install.packages("knitr")
# install.packages("readr")
# install.packages("tikzDevice")
# install.packages("raster")
# install.packages("rgdal")
# install.packages("sf")
# install.packages("rgeos")
# install.packages("gridExtra")
# install.packages("ggvoronoi")
# install.packages("grid")
# install.packages("lattice")
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("scico")
# install.packages("kableExtra")
# install.packages("magrittr")
# install.packages("plyr")
# install.packages("SDraw")

# LOAD PACKAGES
library(Rcpp)
library(WaveSampling)
library(sampling)
library(BalancedSampling)
library(spatstat)
library(spsurvey)
library(sp)
library(SDraw)

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


# knit_hooks$set(document = function(x) {sub('\\usepackage[]{color}', '\\usepackage{xcolor}', x, fixed = TRUE)})
# opts_chunk$set(fig.path='figure/image-', cache.path='cache/latex-')
# options(tikzDocumentDeclaration = "\\documentclass{article}")

pathresults <- file.path(getwd(),"results/",fsep ="/")
if(!dir.exists(pathresults)){
   dir.create(pathresults)
}





## ----GRTS,echo=FALSE---------------------------------------------------------------------------------------------------------------

source("./hippoint2.R")

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

# for hip.point function -> find the index based on coordinates
findIndex <- function(x1,X){
  final <- c()
  for(i in 1:nrow(x1)){
    id <- seq(1,nrow(X),1)
    for(j in 1:ncol(x1)){
      out = which(x1[i,j] == X[id,j])
      if(length(out) > 1){
        id <- out
      }else if(length(out) == 1 & j == 1){
        final <- c(final,out)
        break;
      }else if(length(out) == 1){
        final <- c(final,id[out])
        break;
      }
    }
  }
  return(final)
}



## ----distance,echo = FALSE,warning=FALSE,message = FALSE,results='hide',cache=TRUE-------------------------------------------------

# spatial configuration
X <- expand.grid(seq(0.5,3.5,0.05),seq(0.5,3.5,0.05))
Xp <- expand.grid(1:3,1:3)
N <- nrow(X)

# empty array for the distance
d <- array(rep(0,N*N),c(N,N))
d2 <-  array(rep(0,N*N),c(N,N))
d3 <-  array(rep(0,N*N),c(N,N))

index11 <- 621

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
Xdat <- data.frame(x = X[,1],y = X[,2],v = d[index11,],typ = rep("Tore",N))
Xdat <- rbind(Xdat,data.frame(x = X[,1],y = X[,2],v = d2[index11,],typ = rep("Euclidean",N)))
Xdat <- rbind(Xdat,data.frame(x = X[,1],y = X[,2],v = d3[index11,],typ = rep("Shifted Tore",N)))
Xdat$typ <- factor(Xdat$typ, levels = c("Euclidean", "Tore", "Shifted Tore"))

Xp <- data.frame(x = Xp[,1],y = Xp[,2])
perturb <- data.frame(x = perturb1+1,y = perturb2+1,v = 1,typ = "Shifted Tore")

# ggplot

# tikz(file = "distance.tex", width =   5.78851, height = 3,standAlone = FALSE)
p_distance <- ggplot(Xdat) +
  geom_tile(aes(x = x,y = y,fill = v))+
  geom_point(data = Xp,aes(x = x,y = y),size = 1,pch = 1,stroke = 0.5,colour = "grey50")+
  geom_point(data = perturb,aes(x = x,y = y,fill = v),colour = "black",size = 1)+
  facet_wrap(~typ)+
  coord_fixed()+
  scale_x_continuous(name="$x$",breaks = c(1,2,3) ,limits=c(0.45, 3.55)) +
  scale_y_continuous(name="$y$",breaks = c(1,2,3) ,limits=c(0.45, 3.55)) +
  labs(fill="$m^2$")+
  scale_fill_gradient(low = "white", high = "black",
                     space = "Lab", na.value = "grey50", guide = "colourbar",
                     aesthetics = "fill")+
  theme_wave()
print(p_distance)
# dev.off()


## ----stratification,echo = FALSE,warning=FALSE,message = FALSE,results='hide',cache = TRUE-----------------------------------------

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

# tikz(file = "strat.tex", width =  5.78851, height = 3,standAlone = FALSE)
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


## ----Wsp2,echo = FALSE,warning=FALSE,results='hide',cache = TRUE-------------------------------------------------------------------

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
s_moran <- wave(as.matrix(X_moran),pik_moran,1)

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

# tikz(file = "stratmat.tex", width =  5.78851, height = 3,standAlone = FALSE)
grid.arrange(p_1, p_2, nrow = 1)
# dev.off()



## ----voro1,echo=FALSE, warning=FALSE,message=FALSE,results='hide'------------------------------------------------------------------
	
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

#ggplot

# tikz(file = "voro.tex", width =  5.78851, height = 3,standAlone = FALSE)
p_voron <- ggplot() +
  geom_voronoi(data = vor[vor$s == 1,],aes(x = x,y = y,fill = v),outline =grid,size = 0.1,alpha = 1,colour ="black")+
  geom_point(data = vor,aes(x = x,y = y),size = 1,shape = 1,alpha = 0.5)+
  geom_point(data = vor[vor$s == 1,] ,aes(x = x,y = y),size = 1,shape = 16,alpha = 1)+
  # scale_fill_scico(palette = "vik",begin = 0.7,end = 0.0,limits = c(min(pik),max(sb1)))+
   scale_fill_gradient(low = "grey95", high = "grey20",
                   space = "Lab", na.value = "grey50", guide = "colourbar",
                   aesthetics = "fill")+
  facet_wrap(~name)+
  coord_fixed()+
  labs(fill = "$v_k$")+
  theme_wave()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.border=element_blank())
	print(p_voron)
	# dev.off()


## ----datageneration,echo=FALSE,warning=FALSE,message=FALSE-------------------------------------------------------------------------

	
	# Data loading/generation--------
	
N <- 144
# CSR
# pp <- rpoispp(N)
# X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
# N_pp <- nrow(X_pp)
# while(N_pp != 144){
#   pp <- rpoispp(N)
#   X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
#   N_pp <- nrow(X_pp)
# }
# plot(X_pp)
# save(X_pp,file = "./results/X_pp_20_05_2020.rda")
load("./results/X_pp_20_05_2020.rda")
N_pp <- nrow(X_pp)


# Neyman-Scottt
# nclust <- function(x0, y0, radius, n) {
# return(runifdisc(n, radius, centre=c(x0, y0)))
# }
#
# ns <- rNeymanScott(12, 1, nclust, radius=0.055, n=12)
# X_ns <- matrix(cbind(ns$x,ns$y),nrow = ns$n,ncol = 2)
# N_ns <- nrow(X_ns)
# while(N_ns != 144){
#   print(N_ns)
#   ns <- rNeymanScott(12, 1, nclust, radius=0.055, n=12)
#   X_ns <- matrix(cbind(ns$x,ns$y),nrow = ns$n,ncol = 2)
#   N_ns <- nrow(X_ns)
# }
# plot(X_ns)
# save(X_ns,file = "./results/X_ns_20_05_2020.rda")
load("./results/X_ns_20_05_2020.rda")
N_ns <- nrow(X_ns)




# grid
x <- seq(1,sqrt(N),1)
X_grid <- as.matrix(expand.grid(x,x))
N_grid <- N





## ----artificialPlot,echo = FALSE, warning=FALSE,message=FALSE,results='hide'-------------------------------------------------------

# in order to have the same plot -> could be changed
set.seed(2)

n = 48

# pik
pik_pp <- rep(n/N_pp,N_pp)
pik_ns <- rep(n/N_ns,N_ns)
pik_grid <- rep(n/N,N)

# samples
s_pp <- wave(X_pp,pik_pp,tore = FALSE,shift = FALSE)
s_ns <- wave(X_ns,pik_ns,tore = FALSE,shift = FALSE)
s_grid <- wave(X_grid,pik_grid,tore = TRUE,shift =TRUE)


# data.frame
dat1 <- data.frame(x = X_pp[,1],y = X_pp[,2],type = rep("Complete spatial randomness",nrow(X_pp)))
dat2 <- data.frame(x = X_ns[,1],y = X_ns[,2],type = rep("Neyman-Scott process",nrow(X_ns)))
dat3 <- data.frame(x = X_grid[,1],y = X_grid[,2],type = rep("Simple regular grid",nrow(X_grid)))
dat <- rbind(dat1,dat2,dat3)

# ggplot

# tikz(file = "artificial.tex", width =  5.78851, height = 2.1,standAlone = FALSE)
  p_simul <- ggplot()+
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
  print(p_simul)
# dev.off()


## ----SIMU_pp,echo=FALSE,message=FALSE,results='hide',cache = TRUE,warning=FALSE----------------------------------------------------

# Data loading/generation--------
  
SIM <- 10000
n1 <- 16
n2 <- 32
n3 <- 48
n <- c(n1,n2,n3)

# # pik_pp <- list(rep(n1/N_pp,N_pp),
# #                rep(n2/N_pp,N_pp),
# #                rep(n3/N_pp,N_pp),
# #                inclusionprobabilities(runif(N_pp),n1),
# #                inclusionprobabilities(runif(N_pp),n2),
# #                inclusionprobabilities(runif(N_pp),n3))
# #
# #
# # pik_ns <- list(rep(n1/N_ns,N_ns),
# #                rep(n2/N_ns,N_ns),
# #                rep(n3/N_ns,N_ns),
# #                inclusionprobabilities(runif(N_ns),n1),
# #                inclusionprobabilities(runif(N_ns),n2),
# #                inclusionprobabilities(runif(N_ns),n3))
# #
# # pik_grid <- list(rep(n1/N_grid,N_grid),
# #                  rep(n2/N_grid,N_grid),
# #                  rep(n3/N_grid,N_grid),
# #                  inclusionprobabilities(runif(N_grid),n1),
# #                  inclusionprobabilities(runif(N_grid),n2),
# #                  inclusionprobabilities(runif(N_grid),n3))
# #
# # pik <- list(pik_pp,pik_ns,pik_grid)
# # names(pik) <- c("pp","ns","grid")
# # saveRDS(pik, file = "./results/pik_20_05_2020.rds")
# pik <- readRDS(file = "./results/pik_20_05_2019.rds")
#
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
#     Ib1[[i]][[j]] <- Ib[[i]][[j]] <- Sb[[i]][[j]] <- data.frame(matrix(rep(0,6*SIM),nrow = SIM,ncol =6))
#     colnames(Ib1[[i]][[j]]) <- colnames(Ib[[i]][[j]]) <- colnames(Sb[[i]][[j]]) <-  c("wave","lpm1","scps","grts","hip","srswor")
#
#     s[[i]][[j]] <-  rep(list(0),6)
#     names(s[[i]][[j]]) <- c("wave","lpm1","scps","grts","hip","srswor")
#     # pikt <- UPMEpiktildefrompik(pik[[i]][[j]])
#     # w_maxent <- pikt/(1-pikt)
#     # q[[i]][[j]] <- UPMEqfromw(w_maxent,round(sum(pik[[i]][[j]])))
#
#   }
# }
#
#
# for(tt in 1:SIM){
#   cat("ITERATION ",tt,"\n\n")
#   start_time <- Sys.time()
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
#       # start_time <- Sys.time()
#       s[[i]][[j]][[1]] <- as.vector(wave(X[[i]],pik[[i]][[j]],bound = 1,tore = F,shift = F))
#       # end_time <- Sys.time()
#       # cat("wave time :", end_time - start_time,"\n")
#
#       #lpm1
#       # start_time <- Sys.time()
#       s[[i]][[j]][[2]] <- rep(0,N[[i]])
#       s[[i]][[j]][[2]][lpm1(pik[[i]][[j]],X[[i]])] <- 1
#       # end_time <- Sys.time()
#       # cat("lpm1 time :", end_time - start_time,"\n")
#
#       #scps
#       # start_time <- Sys.time()
#       s[[i]][[j]][[3]] <- rep(0,N[[i]])
#       s[[i]][[j]][[3]][scps(pik[[i]][[j]],X[[i]])] <- 1
#       # end_time <- Sys.time()
#       # cat("scps time :", end_time - start_time,"\n")
#
#       #GRTS
#       # start_time <- Sys.time()
#       s[[i]][[j]][[4]] <- GRTS(pik[[i]][[j]],X[[i]][,1],X[[i]][,2])
#       # end_time <- Sys.time()
#       # cat("GRTS time :", end_time - start_time,"\n")
#
#       #hip
#       if(i != 3 & (j == 1| j == 2| j == 3)){
#         # start_time <- Sys.time()
#         X_Sp <- data.frame(x = X[[i]][,1],y = X[[i]][,2])
#         X_Sp  <- SpatialPoints(X_Sp)
#
#         tmp1 = try(hip.point2(x = X_Sp,
#                               n = sum(pik[[i]][[j]]),
#                               plot.lattice = FALSE),
#                    silent = TRUE)
#         while(class(tmp1)=="try-error"){
#           tmp1 = try(hip.point2(x = X_Sp,
#                              n = sum(pik[[i]][[j]]),
#                              plot.lattice = FALSE),
#                   silent = TRUE)
#         }
#         tmp1 <- coordinates(tmp1)
#         s_tmp <- rep(0,N[[i]])
#         s_tmp[findIndex(tmp1,X[[i]])] <- 1
#         s[[i]][[j]][[5]] <- s_tmp
#         # end_time <- Sys.time()
#         # cat("HIP time :", end_time - start_time,"\n")
#       }else{
#         s[[i]][[j]][[5]] <- 0
#       }
#
#
#       #srswor
#       # start_time <- Sys.time()
#       if(j == 1 | j == 2 | j == 3){
#         s[[i]][[j]][[6]] <- srswor(sum(pik[[i]][[j]]),N[[i]])
#       }else{
#         s[[i]][[j]][[6]] <- UPmaxentropy(pik[[i]][[j]])
#         # s[[j]][[5]] <- UPMEsfromq(q[[j]])
#       }
#       # end_time <- Sys.time()
#       # cat("srswor time :", end_time - start_time,"\n")
#       # cat("\n\n")
#
#     }
#   }
#
#   # fill Ib and Sb
#   for(i in 1:length(s)){
#     for(j in 1:length(s[[i]])){
#       for(k in 1:length(s[[i]][[j]])){
#         if(i == 1 | i == 2){ # pp and ns
#           if(j == 1 | j == 2 | j == 3){ # n1 n2 n3
#             Ib1[[i]][[j]][tt,k] <- IB(W1[[i]][[j]],s[[i]][[j]][[k]])
#             Ib[[i]][[j]][tt,k] <- IB(W[[i]][[j]],s[[i]][[j]][[k]])
#             Sb[[i]][[j]][tt,k] <- sb(pik[[i]][[j]],X[[i]],which(s[[i]][[j]][[k]] == 1))
#           }else{
#             if(k != 5){ # n1_un n2_un n3_un -> no hip
#               Ib1[[i]][[j]][tt,k] <- IB(W1[[i]][[j]],s[[i]][[j]][[k]])
#               Ib[[i]][[j]][tt,k] <- IB(W[[i]][[j]],s[[i]][[j]][[k]])
#               Sb[[i]][[j]][tt,k] <- sb(pik[[i]][[j]],X[[i]],which(s[[i]][[j]][[k]] == 1))
#             }
#           }
#         }else{ #grid
#           if(k != 5){ # no hip
#               Ib1[[i]][[j]][tt,k] <- IB(W1[[i]][[j]],s[[i]][[j]][[k]])
#               Ib[[i]][[j]][tt,k] <- IB(W[[i]][[j]],s[[i]][[j]][[k]])
#               Sb[[i]][[j]][tt,k] <- sb(pik[[i]][[j]],X[[i]],which(s[[i]][[j]][[k]] == 1))
#           }
#         }
#       }
#     }
#   }
#   end_time <- Sys.time()
#   print(end_time - start_time)
#   cat("\n\n")
# }
#
#
# #---- SAVING
#
# pathresults <- "./results"
#
# save <- function(res,name){
#   for(i in 1:length(Ib)){
#     tmp <- do.call(rbind,res[[i]])
#     rownames(tmp) <- NULL
#     tmp$type <- rep(c("n1","n2","n3","n1_un","n2_un","n3_un"),each = SIM)
#     write.table(tmp,
#                 # file = file.path("./results",paste0(name,names(res)[i],".csv")),
#                 file = file.path(pathresults,paste0(name,names(res)[i],".csv")),
#                 row.names=FALSE, na="",col.names=TRUE, sep=",")
#   }
# }
# save(Ib,"Ib")
# save(Ib1,"Ib1")
# save(Sb,"Sb")

# Data loading/generation--------

Ib <- list()

Ibpp <- read_csv("./results/Ibpp_22_05_2020.csv")
Ibpp$type <- factor(Ibpp$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ibpp <- split(as.data.frame(Ibpp),f = Ibpp$type, drop = TRUE)
Ibpp <- lapply(Ibpp, function(x) { x["type"] <- NULL; x })
Ibns <- read_csv("./results/Ibns_22_05_2020.csv")
Ibns$type <- factor(Ibns$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ibns <- split(as.data.frame(Ibns),f = Ibns$type, drop = TRUE)
Ibns <- lapply(Ibns, function(x) { x["type"] <- NULL; x })
Ibgrid <- read_csv("./results/Ibgrid_22_05_2020.csv")
Ibgrid$type <- factor(Ibgrid$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ibgrid <- split(as.data.frame(Ibgrid),f = Ibgrid$type, drop = TRUE)
Ibgrid <- lapply(Ibgrid, function(x) { x["type"] <- NULL; x })

Ib[[1]] <- Ibpp
Ib[[2]] <- Ibns
Ib[[3]] <- Ibgrid

names(Ib) <- c("pp","ns","grid")


Ib1 <- list()

Ib1pp <- read_csv("./results/Ib1pp_22_05_2020.csv")
Ib1pp$type <- factor(Ib1pp$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib1pp <- split(as.data.frame(Ib1pp),f = Ib1pp$type, drop = TRUE)
Ib1pp <- lapply(Ib1pp, function(x) { x["type"] <- NULL; x })
Ib1ns <- read_csv("./results/Ib1ns_22_05_2020.csv")
Ib1ns$type <- factor(Ib1ns$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib1ns <- split(as.data.frame(Ib1ns),f = Ib1ns$type, drop = TRUE)
Ib1ns <- lapply(Ib1ns, function(x) { x["type"] <- NULL; x })
Ib1grid <- read_csv("./results/Ib1grid_22_05_2020.csv")
Ib1grid$type <- factor(Ib1grid$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib1grid <- split(as.data.frame(Ib1grid),f = Ib1grid$type, drop = TRUE)
Ib1grid <- lapply(Ib1grid, function(x) { x["type"] <- NULL; x })

Ib1[[1]] <- Ib1pp
Ib1[[2]] <- Ib1ns
Ib1[[3]] <- Ib1grid

names(Ib1) <- c("pp","ns","grid")

Sb <- list()

Sbpp <- read_csv("./results/Sbpp_22_05_2020.csv")
Sbpp$type <- factor(Sbpp$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Sbpp <- split(as.data.frame(Sbpp),f = Sbpp$type, drop = TRUE)
Sbpp <- lapply(Sbpp, function(x) { x["type"] <- NULL; x })
Sbns <- read_csv("./results/Sbns_22_05_2020.csv")
Sbns$type <- factor(Sbns$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Sbns <- split(as.data.frame(Sbns),f = Sbns$type, drop = TRUE)
Sbns <- lapply(Sbns, function(x) { x["type"] <- NULL; x })
Sbgrid <- read_csv("./results/Sbgrid_22_05_2020.csv")
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
  colnames(Ib1[[i]]) <-colnames(Ib[[i]]) <- colnames(Sb[[i]]) <-  c("design","n","wave","lpm1","scps","grts","hip","srswor","wave","lpm1","scps","grts","hip","maxent")

}




## ----tablepp,echo=FALSE------------------------------------------------------------------------------------------------------------

#data preparation remove hip unequal == 0
datpp <- rbind(Ib1$pp,Ib$pp,Sb$pp)
datpp <- datpp[,-1]
datpp <- datpp[,-12]


colnames(datpp) <- c(" ","wave","lpm1","scps","grts","hip","srswor","wave","lpm1","scps","grts","maxent")

kable(datpp, format = "latex",digits = 3, booktabs = T, caption = paste0("Spreading measures results based on ",SIM," simulations on the Complete spatial randomness dataset. The population size is equal to ",N_pp,"."),row.names = FALSE,escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 6, "Unequal probabilities" = 5),escape = F) %>%
 add_header_above(c(" " = 1,"Sampling design" = 11)) %>%
 group_rows("$I_{B_1}$",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("$I_B$",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("$B$",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9,latex_options="scale_down")



## ----table1,echo=FALSE-------------------------------------------------------------------------------------------------------------

#data preparation remove hip unequal == 0
datIb1 <- rbind(Ib1[[1]],Ib1[[2]],Ib1[[3]])[,2:ncol(Ib1[[1]])]
datIb1 <- datIb1[,-12]
colnames(datIb1) <- c("$I_{B_1}$","wave","lpm1","scps","grts","hip","srswor","wave","lpm1","scps","grts","maxent")

tableIb1 <-  kable(datIb1, format = "latex",digits = 3, booktabs = T, caption = paste0("Results of ",SIM," simualtions of Moran's $I_{B_1}$ spatial measure \\eqref{eq:wpik1}.","The three spatial configurations of the Section \\ref{sec:artificial} are taken into account."),escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 6, "Unequal probabilities" = 5),escape = F) %>%
 add_header_above(c(" " = 1,"Sampling design" = 11)) %>%
 group_rows("Complete spatial randomness",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Neyman-Scott process",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Simple regular grid",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9,latex_options="scale_down")


## ----table2,echo=FALSE-------------------------------------------------------------------------------------------------------------

#data preparation remove hip unequal == 0
datIb <- rbind(Ib[[1]],Ib[[2]],Ib[[3]])[,2:ncol(Ib[[1]])]
datIb <- datIb[,-12]
colnames(datIb) <- c("$I_{B}$","wave","lpm1","scps","grts","hip","srswor","wave","lpm1","scps","grts","maxent")

tableIb <- kable(datIb, format = "latex",digits = 3, booktabs = T, caption = paste0("Results of ",SIM," simulations of Moran's $I_B$ spatial measure \\eqref{eq:wpik}.","The three spatial configurations of the Section \\ref{sec:artificial} are taken into account."),row.names = FALSE,escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 6, "Unequal probabilities" = 5)) %>%
 add_header_above(c(" " = 1,"Sampling design" = 11)) %>%
 group_rows("Complete spatial randomness",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Neyman-Scott process",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Simple regular grid",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9,latex_options="scale_down")



## ----table3,echo=FALSE-------------------------------------------------------------------------------------------------------------

#data preparation remove hip unequal == 0

datSb <- rbind(Sb[[1]],Sb[[2]],Sb[[3]])[,2:ncol(Sb[[1]])]
datSb <- datSb[,-12]
colnames(datSb) <- c("$B$","wave","lpm1","scps","grts","hip","srswor","wave","lpm1","scps","grts","maxent")

tableSb <- kable(datSb, format = "latex",digits = 3, booktabs = T, caption = paste0("Results of ",SIM," simulations of the spatial balance measure $B$ \\eqref{eq:voro}.","The three spatial configurations of the Section \\ref{sec:artificial} are taken into account."),row.names = FALSE,escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 6, "Unequal probabilities" = 5)) %>%
 add_header_above(c(" " = 1,"Sampling design" = 11)) %>%
 group_rows("Complete spatial randomness",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Neyman-Scott process",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("Simple regular grid",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9,latex_options="scale_down")




## ----meuse_plot,echo =FALSE--------------------------------------------------------------------------------------------------------
data("meuse")
data("meuse.riv")
meuse.riv <- meuse.riv[which(meuse.riv[,2] < 334200 & meuse.riv[,2] > 329400),]

X <- scale(as.matrix(meuse[,1:2]))
X_SpatialPoint <- data.frame(x = X[,1],y = X[,2])
X_SpatialPoint  <- SpatialPoints(X_SpatialPoint, proj4string=CRS("+init=epsg:28992"), bbox = NULL)

pik <- inclusionprobabilities(meuse$copper,30)


## ----meusefig,echo = FALSE, warning=FALSE,message=FALSE,results='hide',cache = TRUE------------------------------------------------

meuse_sf <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
s <- wave(X,pik)

# tikz(file = "meuse.tex", width =5, height = 4,standAlone = FALSE)
meuseplot <- ggplot()+
  geom_sf(data = meuse_sf,aes(size=copper),show.legend = 'point',shape = 1,stroke = 0.3)+
  geom_polygon(data = data.frame(x = meuse.riv[,1],y = meuse.riv[,2]),
               aes(x = x,y = y),
               fill = "grey70",
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



## ----spread_measure_meuse,echo = FALSE,message=FALSE,warning=FALSE,results='hide',cache = TRUE-------------------------------------


SIM <- 10000
n1 <- 15
n2 <- 30
n3 <- 50
n <- c(n1,n2,n3)
N_meuse <- nrow(X)
y <- meuse$cadmium*0.01


pik_meuse <- list(rep(n1/N_meuse,N_meuse),
                  rep(n2/N_meuse,N_meuse),
                  rep(n3/N_meuse,N_meuse),
                  inclusionprobabilities(meuse$copper,n1),
                  inclusionprobabilities(meuse$copper,n2),
                  inclusionprobabilities(meuse$copper,n3))


# #---- INITIALIZATION EMPTY CONTAINER
#
# s <- W1<- W <- Ib <- Ib1 <- Sb  <- y_hat <- sim <- V <- V_ht <- rep(list(0),6)
# names(s) <- names(W) <-names(W1) <- names(Ib1) <- names(Ib) <- names(Sb) <-names(y_hat)<- names(sim) <-names(V) <- names(V_ht) <- c("n1","n2","n3","n1_un","n2_un","n3_un")
# for(j in 1:length(s)){
#   # V_ht[[j]] <-  data.frame(matrix(rep(0,5*SIM),nrow = SIM,ncol = 5))
#   V_ht[[j]] <- rep(list(0),6)
#   names(V_ht[[j]]) <- c("wave","lpm1","scps","grts","hip","srswor")
#   for(k in 1:length(V_ht[[j]])){
#     V_ht[[j]][[k]] <-  data.frame(matrix(rep(0,6*SIM),nrow = SIM,ncol = 6))
#     colnames(V_ht[[j]][[k]]) <- c("VSIM","VSB","VNBH2","VNBH3","VNBH4","VHAJ")
#   }
#
#
#   Ib1[[j]] <- Ib[[j]] <- Sb[[j]] <- sim[[j]] <- V[[j]] <- y_hat[[j]]  <- data.frame(matrix(rep(0,6*SIM),nrow = SIM,ncol = 6))
#   colnames(Ib1[[j]]) <- colnames(Ib[[j]]) <- colnames(Sb[[j]]) <- colnames(sim[[j]]) <- colnames(V[[j]]) <- colnames(y_hat[[j]]) <-  c("wave","lpm1","scps","grts","hip","srswor")
#   s[[j]] <-  rep(list(0),6)
#   names(s[[j]]) <- c("wave","lpm1","scps","grts","hip","srswor")
# }
#
#
# countTest <- 1
# #---- LOOP SIMULATION
# for(tt in 1:SIM){
#   cat("\n ITERATION ",tt,"\n\n")
#   for(j in 1:length(s)){
#     #W
#     W[[j]] <- wpik(X,pik_meuse[[j]],tore = F,shift = F)
#     W[[j]]  <- drop0(W[[j]] - diag(diag(W[[j]])), tol = 0, is.Csparse = NA)
#     #W1
#     W1[[j]] <- wpikInv(X,pik_meuse[[j]],tore = F,shift = F,toreBound = 0)
#
#     #wave
#     # start_time <- Sys.time()
#     s[[j]][[1]] <- as.vector(wave(X,pik_meuse[[j]],bound = 1,tore = F,shift = F))
#     # end_time <- Sys.time()
#     # cat("WAVE time :", end_time - start_time,"\n")
#
#     #lpm1
#     # start_time <- Sys.time()
#     s[[j]][[2]] <- rep(0,N_meuse)
#     s[[j]][[2]][lpm1(pik_meuse[[j]],X)] <- 1
#     # end_time <- Sys.time()
#     # cat("LPM time :", end_time - start_time,"\n")
#
#
#     #scps
#     # start_time <- Sys.time()
#     s[[j]][[3]] <- rep(0,N_meuse)
#     s[[j]][[3]][scps(pik_meuse[[j]],X)] <- 1
#     # end_time <- Sys.time()
#     # cat("SCPS time :", end_time - start_time,"\n")
#
#     #grts
#     # start_time <- Sys.time()
#     s[[j]][[4]] <- GRTS(pik_meuse[[j]],X[,1],X[,2])
#     # end_time <- Sys.time()
#     # cat("GRTS time :", end_time - start_time,"\n")
#
#     #hip
#     if(j == 1| j == 2| j == 3){
#       # start_time <- Sys.time()
#
#       tmp1 = try(hip.point2(x = X_SpatialPoint,
#                             n = sum(pik_meuse[[j]]),
#                             plot.lattice = FALSE),
#                  silent = TRUE)
#       while(class(tmp1)=="try-error"){
#         tmp1 = try(hip.point2(x = X_SpatialPoint,
#                            n = sum(pik_meuse[[j]]),
#                            plot.lattice = FALSE),
#                 silent = TRUE)
#         countTest <- countTest + 1
#       }
#       tmp1 <- coordinates(tmp1)
#       s_tmp <- rep(0,N_meuse)
#       s_tmp[findIndex(tmp1,X)] <- 1
#       s[[j]][[5]] <- s_tmp
#       # end_time <- Sys.time()
#       # cat("HIP time :", end_time - start_time,"\n")
#     }else{
#       s[[j]][[5]] <- 0
#     }
#
#
#
#     #srswor
#     # start_time <- Sys.time()
#     if(j == 1 | j == 2 | j == 3){
#       s[[j]][[6]] <- srswor(sum(pik_meuse[[j]]),N_meuse)
#     }else{
#       s[[j]][[6]] <- UPmaxentropy(pik_meuse[[j]])
#       # s[[j]][[5]] <- UPMEsfromq(q[[j]])
#     }
#     # end_time <- Sys.time()
#     # cat("SRS (MAXENT) time :", end_time - start_time,"\n")
#     # cat("\n\n")
#
#
#   }
#
#   # fill Ib and Sb
#   for(j in 1:length(s)){
#     for(k in 1:length(s[[j]])){
#       if((j == 4 | j == 5 | j == 6) & k == 5){
#
#       }else{
#         Ib1[[j]][tt,k] <- IB(W1[[j]],s[[j]][[k]])
#         Ib[[j]][tt,k] <- IB(W[[j]],s[[j]][[k]])
#         Sb[[j]][tt,k] <- sb(pik_meuse[[j]],X,which(s[[j]][[k]] == 1))
#       }
#     }
#   }
#
#
#   # variance
#   for(j in 1:length(s)){
#     for(k in 1:length(s[[j]])){
#       if((j == 4 | j == 5 | j == 6) & k == 5){
#
#       }else{
#         y_hat[[j]][tt,k] <- sum(y[which(s[[j]][[k]]==1)]/pik_meuse[[j]][which(s[[j]][[k]] ==1)])
#         sim[[j]][tt,k] <- (sum(y[which(s[[j]][[k]]==1)]/pik_meuse[[j]][which(s[[j]][[k]] ==1)]) - sum(y))^2
#
#         if(names(s[[j]])[k] == "srswor"){
#           V[[j]][tt,k] <- varHAJ(y[which(s[[j]][[k]]==1)],pik_meuse[[j]][which(s[[j]][[k]]==1)],which(s[[j]][[k]]==1))
#         }else if(names(s[[j]])[k] == "wave"){
#           V[[j]][tt,k] <- BalancedSampling::vsb(probs = pik_meuse[[j]][which(s[[j]][[k]]==1)],
#                                                 ys = y[which(s[[j]][[k]]==1)],
#                                                 xs = X[which(s[[j]][[k]]==1),])
#           # localmean.var(y[which(s[[j]][[k]]==1)]/pik_meuse[[j]][which(s[[j]][[k]]==1)],wtmp)
#
#           # wtmp <- localmean.weight(x = X[which(s[[j]][[k]]==1),1],
#           #                          y = X[which(s[[j]][[k]]==1),2],
#           #                          prb = pik_meuse[[j]][which(s[[j]][[k]]==1)],
#           #                          nbh = 2)
#           # V[[j]][tt,k] <- localmean.var(y[which(s[[j]][[k]]==1)]/pik_meuse[[j]][which(s[[j]][[k]]==1)],wtmp)
#         }else{
#           V[[j]][tt,k] <- BalancedSampling::vsb(probs = pik_meuse[[j]][which(s[[j]][[k]]==1)],
#                                                 ys = y[which(s[[j]][[k]]==1)],
#                                                 xs = X[which(s[[j]][[k]]==1),])
#           # wtmp <- localmean.weight(x = X[which(s[[j]][[k]]==1),1],
#           #                          y = X[which(s[[j]][[k]]==1),2],
#           #                          prb = pik_meuse[[j]][which(s[[j]][[k]]==1)],
#           #                          nbh = 4)
#           # V[[j]][tt,k] <- localmean.var(y[which(s[[j]][[k]]==1)]/pik_meuse[[j]][which(s[[j]][[k]]==1)],wtmp)
#         }
#       }
#     }
#   }
#
#
#   # different variance
#   for(j in 1:length(s)){
#     for(k in 1:length(V_ht[[j]])){
#
#       if((j == 4 | j == 5 | j == 6) & k == 5){
#
#       }else{
#         # vsim
#         V_ht[[j]][[k]][tt,1] <- (sum(y[which(s[[j]][[k]] == 1)]/pik_meuse[[j]][which(s[[j]][[k]] ==1)]) - sum(y))^2
#
#         # VSB
#         V_ht[[j]][[k]][tt,2] <- BalancedSampling::vsb(probs = pik_meuse[[j]][which(s[[j]][[k]] == 1)],
#                                                  ys = y[which(s[[j]][[k]] == 1)],
#                                                  xs = X[which(s[[j]][[k]] == 1),])
#         # VNBH2
#         wtmp <- localmean.weight(x = X[which(s[[j]][[k]] == 1),1],
#                                  y = X[which(s[[j]][[k]] == 1),2],
#                                  prb = pik_meuse[[j]][which(s[[j]][[k]] == 1)],
#                                  nbh = 2)
#         V_ht[[j]][[k]][tt,3] <- localmean.var(y[which(s[[j]][[k]] == 1)]/pik_meuse[[j]][which(s[[j]][[k]] == 1)],wtmp)
#         # VNBH3
#         wtmp <- localmean.weight(x = X[which(s[[j]][[k]] == 1),1],
#                                  y = X[which(s[[j]][[k]] == 1),2],
#                                  prb = pik_meuse[[j]][which(s[[j]][[k]] == 1)],
#                                  nbh = 3)
#         V_ht[[j]][[k]][tt,4] <- localmean.var(y[which(s[[j]][[k]] == 1)]/pik_meuse[[j]][which(s[[j]][[k]] == 1)],wtmp)
#
#         #VNBH4
#         wtmp <- localmean.weight(x = X[which(s[[j]][[k]] == 1),1],
#                                  y = X[which(s[[j]][[k]] == 1),2],
#                                  prb = pik_meuse[[j]][which(s[[j]][[k]] == 1)],
#                                  nbh = 4)
#         V_ht[[j]][[k]][tt,5] <- localmean.var(y[which(s[[j]][[k]] == 1)]/pik_meuse[[j]][which(s[[j]][[k]] == 1)],wtmp)
#
#         #VHAJ
#         V_ht[[j]][[k]][tt,6] <- varHAJ(y[which(s[[j]][[k]]==1)],pik_meuse[[j]][which(s[[j]][[k]]==1)],which(s[[j]][[k]]==1))
#
#       }
#     }
#   }
#
#
# }
#
#
# #------ V_ht
#
# for(i in 1:length(V_ht)){
#   V_ht[[i]] <- lapply(V_ht[[i]],colMeans)
#   tmp <- c()
#   for(j in 1:length(V_ht[[i]])){
#     tmp <- cbind(tmp,V_ht[[i]][[j]])
#   }
#   colnames(tmp) <- c("wave","lpm1","scps","grts","hip","srswor")
#   V_ht[[i]] <- tmp
# }
#
# V_ht <- do.call(rbind,V_ht)
# V_ht <- as.data.frame(V_ht)
# rownames(V_ht) <- paste(rownames(V_ht),rep(c("n1","n2","n3","n1_un","n2_un","n3_un"),each = 6),sep = ":")
# V_ht <- cbind(V_ht[1:18,],V_ht[19:36,])
# V_ht <- cbind(rep(c("$v_{SIM}$","$v_{SB}$","$v_{LM2}$","$v_{LM3}$","$v_{LM4}$","$v_{HAJ}$"),3),V_ht)
# colnames(V_ht)[1] <- " "
#
#
# #---- SAVING
#
#
# pathresults <- "./results"
#
#
# write.table(V_ht,
#               file =  file.path(pathresults,paste0("V_htMeuse",".csv")),
#               row.names=FALSE, na="",col.names=TRUE, sep=",")
#
#
#
# save <- function(res,name){
#   tmp <- do.call(rbind,res)
#   rownames(tmp) <- NULL
#   tmp$type <- rep(c("n1","n2","n3","n1_un","n2_un","n3_un"),each = SIM)
#   write.table(tmp,
#               file =  file.path(pathresults,paste0(name,".csv")),
#               row.names=FALSE, na="",col.names=TRUE, sep=",")
# }
#
# save(Ib,"IbMeuse")
# save(Ib1,"Ib1Meuse")
# save(Sb,"SbMeuse")
# save(V,"varMeuse")
# save(sim,"VsimuMeuse")
# save(y_hat,"YhtMeuse")




#---- LOADING

V_ht <- read.csv("./results/V_htMeuse_22_05_2020.csv")


Ib <- read_csv("./results/IbMeuse_22_05_2020.csv")
Ib$type <- factor(Ib$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib <- split(as.data.frame(Ib),f = Ib$type, drop = TRUE)
Ib <- lapply(Ib, function(x) { x["type"] <- NULL; x })

Ib1 <- read_csv("./results/Ib1Meuse_22_05_2020.csv")
Ib1$type <- factor(Ib1$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Ib1 <- split(as.data.frame(Ib1),f = Ib1$type, drop = TRUE)
Ib1 <- lapply(Ib1, function(x) { x["type"] <- NULL; x })

Sb <- read_csv("./results/SbMeuse_22_05_2020.csv")
Sb$type <- factor(Sb$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
Sb <- split(as.data.frame(Sb),f = Sb$type, drop = TRUE)
Sb <- lapply(Sb, function(x) { x["type"] <- NULL; x })

V <- read_csv("./results/varMeuse_22_05_2020.csv")
V$type <- factor(V$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
V <- split(as.data.frame(V),f = V$type, drop = TRUE)
V <- lapply(V, function(x) { x["type"] <- NULL; x })

sim <- read_csv("./results/VsimuMeuse_22_05_2020.csv")
sim$type <- factor(sim$type,levels = c("n1","n2","n3","n1_un","n2_un","n3_un"))
sim <- split(as.data.frame(sim),f = sim$type, drop = TRUE)
sim <- lapply(sim, function(x) { x["type"] <- NULL; x })

y_hat <- read_csv("./results/YhtMeuse_22_05_2020.csv")
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

Ib1 <- Ib1[,-12]
Ib  <- Ib[,-12]
Sb  <- Sb[,-12]

colnames(Ib1) <-colnames(Ib) <- colnames(Sb) <-  c(" ","wave","lpm1","scps","grts","hip","srswor","wave","lpm1","scps","grts","maxent")





## ----spread_meuse_table,echo=FALSE,results='asis'----------------------------------------------------------------------------------

datmeuse <- rbind(Ib1,Ib,Sb)
# saveRDS(datvar, file = "C:/Users/jauslinr/switchdrive/SpreadSampling/draft/table/SPR.rds")


kable(datmeuse, format = "latex",digits = 3, booktabs = T, caption = paste0("Spreading measures results based on ",SIM," simulations on the Meuse dataset. The population size is equal to ",N_meuse,"."),row.names = FALSE,escape = FALSE)%>%
 add_header_above(c(" " = 1,"Equal probabilities" = 6, "Unequal probabilities" = 5),escape = F) %>%
 add_header_above(c(" " = 1,"Sampling design" = 11)) %>%
 group_rows("$I_{B_1}$",1,3,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("$I_B$",4,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
 group_rows("$B$",7,9,escape = F,bold = F,latex_gap_space = "1ex")%>%
  kable_styling(font_size = 9,latex_options="scale_down")





## ----var_meuse,echo = FALSE,message=FALSE,warning=FALSE,results='hide'-------------------------------------------------------------

# Coverage rate
Up <- Low <- y_hat
for(l in 1:6){
  for(r in 1:6){
    Up[[l]][,r] <- y_hat[[l]][,r] + qnorm(0.5 + (95/100)/2)*sqrt(V[[l]][,r])
    Low[[l]][,r] <- y_hat[[l]][,r] - qnorm(0.5 + (95/100)/2)*sqrt(V[[l]][,r])
  }
}

coverage <- Up
for(l in 1:6){
  for(r in 1:6){
    coverage[[l]][,r] <- (Low[[l]][,r] < sum(y) &  sum(y) < Up[[l]][,r])
  }
}

coverage_m <- lapply(coverage,colMeans)
sim_m <- lapply(sim,colMeans)
V_m <- lapply(V,colMeans)

# Ration between averages of estimated and simu
aver_var <- sim_m
for(l in 1:6){
  for(r in 1:6){
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


sim <- sim[,-12]
V <- V[,-12]
coverage_m <- coverage_m[,-12]
aver_var <- aver_var[,-12]


colnames(sim) <-colnames(V) <- colnames(coverage_m) <-colnames(aver_var) <-  c(" ","wave","lpm1","scps","grts","hip","srswor","wave","lpm1","scps","grts","maxent")




## ----tablevar,echo=FALSE-----------------------------------------------------------------------------------------------------------


datvar <- rbind(sim,V,coverage_m,aver_var)
# saveRDS(datvar, file = "../draft/table/VAR.rds")


kable(datvar, format = "latex",digits = 3, booktabs = T, caption = paste0("Results of ",SIM," simulations on Meuse dataset. The population size is equal to ",N_meuse,"."," $v_{SIM}$ is equal to the variance approximated by the simulations \\eqref{eq:varSIM}. $v$ depends on the sampling design.
For the srswor and maxent methods, we used the estimator $v_{HAJ}$ \\eqref{eq:varHAJ} while for the other sampling designs, we use $v_{SB}$ \\eqref{eq:varSB}. Coverage rate of the 95\\% confidence intervals are computed as well as the ratio between averages of $v$ and $v_{SIM}$."),row.names = FALSE,escape = FALSE)%>%
  add_header_above(c(" " = 1,"Equal probabilities" = 6, "Unequal probabilities" = 5),escape = F) %>%
  add_header_above(c(" " = 1,"Sampling design" = 11)) %>%
  group_rows("$v_{SIM}$",1,3,escape = F,bold = F,latex_gap_space = "1em")%>%
  group_rows("$v$",4,6,escape = F,bold = F,latex_gap_space = "1em")%>%
  group_rows("Coverage of the 95% confidence interval",7,9,escape = T,bold = F,latex_gap_space = "1em") %>%
  group_rows("Ratio $v/v_{SIM}$",10,12,escape = F,bold = F,latex_gap_space = "1em") %>%
  kable_styling(font_size = 9,latex_options="scale_down")



## ----V_ht_meuse,echo = FALSE,message=FALSE,warning=FALSE,results='asis'------------------------------------------------------------

V_ht <- V_ht[,-12]
colnames(V_ht) <- c(" ","wave","lpm1","scps","grts","hip","srswor","wave","lpm1","scps","grts","maxent")



kable(V_ht, format = "latex",digits = 3, booktabs = T,
      caption = paste0("Results of ",SIM," simulations on Meuse dataset. The population size is equal to ",N_meuse,"."," $v_{SIM}$ \\eqref{eq:varSIM} is equal to the variance approximated by the simulations. $v_{SB}$ \\eqref{eq:varSB} is the variance estimator based on the nearest neighbours in the sample. $v_{LM_j}$ is equal to the estimator \\eqref{eq:varLM} where the number of neighbouring units used is set to $j = 2,3,4$. $v_{HAJ}$ \\eqref{eq:varHAJ} is the Hajek-Rosen estimator."),
row.names = FALSE,escape = FALSE)%>%
  add_header_above(c(" " = 1,"Equal probabilities" = 6, "Unequal probabilities" = 5),escape = F) %>%
  add_header_above(c(" " = 1,"Sampling design" = 11)) %>%
  group_rows("$n = 15$",1,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
  group_rows("$n = 30$",7,12,escape = F,bold = F,latex_gap_space = "1ex")%>%
  group_rows("$n = 50$",13,18,escape = F,bold = F,latex_gap_space = "1ex")%>% kable_styling(font_size = 9,latex_options="scale_down")



