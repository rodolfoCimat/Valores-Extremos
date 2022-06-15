rm(list =ls())
library(dplyr)
require(ismev); require(SpatialExtremes); require(latex2exp);require(fitdistrplus)
require(reshape2)

##########  Cool Plots  ###########
plotR<-function(x,y=TRUE,GRID = 1,...){
if(class(y) != "logical"){
plot(x,y,...)
if(GRID == 1){
grid()
}
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2) 
}
else{
plot(x,...)
if(GRID == 1){
grid()
}
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2) 
}
}

FEMa<-function(u,X)ifelse(sum(X > u) >0,mean(X[X > u] - u),0)
FEM<-function(u,X)sapply(u,function(i)FEMa(i,X))

ED_estimator <- function(X, k){
  n <- length(X)
  res <- X[n]
  for (i in 1:k){
   a_ik <- log( (k+i)/(k+i - 1.0) )/log(2.0)
   res <- res + a_ik*( X[n - k] - X[n - k - i - 1] )
  }
  return(res)
}

G_nk0 <- function(X, k){
  n <- length(X)
  x_F <- ED_estimator(X,k)
  G_nk <- (x_F - X[n - k])/(X[n-k] - X[n-2*k])
  G <- log(2.0)*G_nk - (log(k) + 0.5*log(2.0))
  
  return(G)
}

test_zero_vs_positive <- function(X,k){
  G <- G_nk0(X,k)
  p <- 1 - exp(-exp(-G))
  return(p)
}

test_zero_vs_negative <- function(X,k){
  G <- G_nk0(X,k)
  p <- exp(-exp(-G))
  return(p)
}

####################################
SeaLev<-read.csv(file.choose()) 
SeaLev<-as_tibble(SeaLev)
SeaLev<-dplyr::select(SeaLev,c("Year","GMSL_noGIA","StdDevGMSL_noGIA"))
names(SeaLev)<-c("year","gmsl","STD")
SeaLev<-mutate(SeaLev, uno = rep(1, nrow(SeaLev)))
head(SeaLev,3)
tail(SeaLev,3)


#######Inciso a) 
par(mfrow = c(1,2))
j<- "gmsl"
(gmsl<- unname(unlist(dplyr::select(SeaLev,all_of(j)))))
gm<-gmsl
plotR(gmsl, type = "p", pch = 20, col = "red", 
       main = paste("Gráfico de dispersión: ",j))
summary(gmsl)
j<- "STD"
(STD<- unname(unlist(dplyr::select(SeaLev,all_of(j)))))
std1<-STD
plotR(STD, type = "p", pch = 20, col = "red", 
       main = paste("Gráfico de dispersión: ",j))


####### Inciso b) Determinemos si existe un dominio de atracción máximal 
#######
###FEP's 
par(mfrow = c(1,1))
gmsl <-sort(gmsl)
ff1<-FEM(gmsl, gmsl)
plotR(gmsl,ff1, type = "l", col = "purple", lwd = 2,
       main = "FEP GMSL")

par(mfrow = c(1,1))
STD <-sort(STD)
ff1<-FEM(STD,STD)
plotR(STD,ff1, type = "l", col = "purple", lwd = 2,
      , main = "FEP STD")###0.010
abline(v = 93)##0.084
abline(v = 100)##0.037
sum(STD > 93);sum(STD > 100);sum(STD > 110)
mean(STD > 93);mean(STD > 100);mean(STD > 110)


###Cociente GMSL No hay evidencia contra Weibull at this points
par(mfrow = c(1,2))
gmsl <-sort(gmsl)
ff1<-FEM(gmsl, gmsl)
plotR(gmsl,ff1/gmsl, type = "l", col = "purple", lwd = 2,
      main = "Cociente FEP empír. GMSL")
u<-quantile(gmsl,0.75)
q<-quantile(gmsl,0.985)
mean(gmsl > u);mean(gmsl > q)
xlim<-c(u,q) 
plotR(gmsl,ff1/gmsl, type = "l", col = "purple", lwd = 2
      ,xlim=xlim,ylim = c(0,0.7),  main = "Cociente FEP empír. GMSL Zoom")

###Cociente FEP STD No evidencia en contra Weibull hasta este punto
STD <-sort(STD)
ff2<-FEM(STD, STD)
plotR(STD,ff2/STD, type = "l", col = "purple", lwd = 2,
      main = "Cociente FEP empír. STD")
(u<-quantile(STD,0.75))
(q<-quantile(STD,0.985))
mean(STD > u);mean(STD > q)
xlim<-c(u,q) 
plotR(STD,ff2/STD,type = "l", col = "purple", lwd = 2
      ,xlim=xlim,main = "Cociente FEP empír STD Zoom")

######################### Datos transformados ######################: 

###FEP transformados
gmsl <-sort(gmsl)
MM<-max(gmsl)
tgmsl<- (MM- gmsl[gmsl != MM])^(-1)
fft1<-FEM(tgmsl,tgmsl)
plotR(tgmsl,fft1, type = "l", col = "purple", lwd = 2,
      main = "FEP empír. GMSL Trans")
u<-quantile(tgmsl,0.65)
q<-quantile(tgmsl,0.985)
mean(tgmsl > u);mean(tgmsl > q)
xlim<-c(u,q) 
plotR(tgmsl,fft1, type = "l", col = "purple", lwd = 2
      ,xlim=xlim,  main = "FEP empír. GMSL Trans Zoom", ylim = c(0,.8))
)

###Cociente FEP transformados
plotR(tgmsl,fft1/tgmsl, type = "l", col = "purple", lwd = 2,
      main = "Cociente FEP empír. GMSL Trans")
u<-quantile(tgmsl,0.65)
q<-quantile(tgmsl,0.985)
mean(tgmsl > u);mean(tgmsl > q)
xlim<-c(u,q) 
plotR(tgmsl,fft1/tgmsl, type = "l", col = "purple", lwd = 2
      ,xlim=xlim,  main = "Cociente FEP empír. GMSL Trans Zoom")

####
###FEP transformados STD
STD <-sort(STD)
MM<-max(STD)
tgmsl<- (MM- STD[STD != MM])^(-1)
fft1<-FEM(tgmsl,tgmsl)
plotR(tgmsl,fft1, type = "l", col = "purple", lwd = 2,
      main = "Cociente FEP empír. STD Trans")
u<-quantile(tgmsl,0.65)
q<-quantile(tgmsl,0.985)
mean(tgmsl > u);mean(tgmsl > q)
xlim<-c(u,q) 
plotR(tgmsl,fft1, type = "l", col = "purple", lwd = 2
      ,xlim=xlim,  main = "FEP empír. STD Trans Zoom", ylim = c(0,0.2))
)

###Cociente FEP transformados STD
plotR(tgmsl,fft1/tgmsl, type = "l", col = "purple", lwd = 2,
      main = "Cociente FEP empír. STD Trans")
u<-quantile(tgmsl,0.65)
q<-quantile(tgmsl,0.985)
mean(tgmsl > u);mean(tgmsl > q)
xlim<-c(u,q) 
plotR(tgmsl,fft1/tgmsl, type = "l", col = "purple", lwd = 2
      ,xlim=xlim,  main = "Cociente FEP empír. STD Trans Zoom")
)

#######################################################################
######################################################Ajustes DGP :DD

(u<-quantile(gmsl,0.94))
(q<-quantile(gmsl,0.95))
mean(gmsl > u);mean(gmsl > q)
gpd.fitrange(gmsl,u,q) 
###Zoom elección umbraaal 
gpd.fitrange(gmsl,48.3,48.6)
gpd.fitrange(gmsl,47.5,47.7) 
u.gmsl<-47.55
mean(gmsl > u.gmsl)
fit.gmsl<-gpd.fit(gmsl,u.gmsl)
gpd.diag(fit.gmsl)
round(fit.gmsl$mle,3) 

summary(STD)
(u<-quantile(STD,0.94))
(q<-quantile(STD,0.95))
mean(STD > u);mean(STD > q)
gpd.fitrange(STD,u,q) 
###Zoom elección umbraaal 
gpd.fitrange(STD,u,97) 
u.STD<-96.5
mean(STD > u.STD)
fit.STD<-gpd.fit(STD,u.STD)
gpd.diag(fit.STD)
round(fit.STD$mle,3) 

#################################################################

#####################################                          Ajuste DGVE
SL <-group_by(SeaLev,year) 
(n <-summarise(SL, uno = sum(uno)))

####Max por bloques de 8.

maxgmsl<- sapply(1:131,function(j)max(gm[(8*(j-1) +1):(8*j)])) 
fit.gmsl <- gev.fit(maxgmsl)
round(fit.gmsl$mle,3)
gev.diag(fit.gmsl)
length(maxgmsl)

####Max por bloques de 8.
maxSTD <- sapply(1:131,function(j)max(std1[(8*(j-1) +1):(8*j)])) 
fit.MSTD <- gev.fit(maxSTD)
round(fit.MSTDl$mle,3)
gev.diag(fit.MSTD)


#####STD G_{n.k}(0)
mon<-sort(STD)
n<-length(mon)
p_negative <- test_zero_vs_negative(mon,1)
      for (k in 2:(floor(n^0.6))){
  		p_negative <- c(p_negative, test_zero_vs_negative(mon,k))
	}

p_positive <- test_zero_vs_positive(mon,1)
      for (k in 2:(floor(n^0.6))){
  		p_positive <- c(p_positive, test_zero_vs_positive(mon,k))
	}

plotR(p_negative,type = "l",lwd = 2, col ="purple",
      main = TeX("STD $\\xi < 0$ vs $\\xi = 0$"))
abline(h = 0.05, col = "red")
plotR(p_positive,type = "l",lwd = 2, col ="purple",
           main = TeX("STD $\\xi > 0$ vs $\\xi = 0$"))
abline(h = 0.05, col = "red")



######################################Estimaciones ED metodología ayudantía

par(mfrow = c(1,2))
mon<-sort(gmsl)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 1:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho GMSL"),  lwd = 2,
     xlab = "n")
(w.gmsl<-mean(w_F))
abline(h =w.gmsl,col= "red")


####
mon<-sort(STD)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 1:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho STD"),  lwd = 2,
     xlab = "n")
(w.STD<-mean(w_F))
abline(h =w.STD,col= "red")
max(STD)

##############################################Periodos de retorno 
###Informal 
xx<- 50 
(colapa <-pgpd(xx-u.gmsl,loc = 0, 
           scale = fit.gmsl$mle[1],shape =fit.gmsl$mle[2], lower.tail = FALSE))
(colaemp <-fit.gmsl$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
(ret.gmsl1<-1/estim)

###Poco más formal 
mon<-gmsl

(w.gmsl<-(w.gmsl + max(mon))/2 )
Mon <-1/(w.gmsl- mon)
Mon <-sort(Mon)
(uu <-1/(w.gmsl- u.gmsl))
(xx <-1/(w.gmsl - xx))
fitgm <-gpd.fit(Mon,uu)
round(fitgm$mle,3)
gpd.diag(fitgm)

(colapa <-pgpd(xx-uu,loc = 0, 
           scale = fitgm$mle[1],shape =fitgm$mle[2], lower.tail = FALSE))
(colaemp <-fitgm$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
(ret.gmsl2<-1/estim)
c(ret.gmsl2,ret.gmsl1)

##Con DGVE 

n <-summarise(SL, uno = sum(uno))
(n<- dplyr::select(n,uno))
summary(n)
(n <- floor(unlist(unname(apply(n, 2,mean)))))
(mle<- fit.Mgmsl$mle)
(dist<-pgev(50,loc = mle[1],scale=mle[2],shape= mle[3]))
1-dist
(estim <- 1 - dist^(1/n))
(ret.gmsl3<-1/estim)


###################################

xx<- 100
(colapa <-pgpd(xx-u.STD,loc = 0, 
           scale = fit.STD$mle[1],shape =fit.STD$mle[2], lower.tail = FALSE))
(colaemp <-fit.STD$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
(ret.STD1<-1/estim)

###Poco más formal 
mon<-STD
w.STD
Mon <-1/(w.STD- mon)
Mon <-sort(Mon)
(uu <-1/(w.STD- u.STD))
(xx <-1/(w.STD - xx))
fitgm <-gpd.fit(Mon,uu)
round(fitgm$mle,3)
gpd.diag(fitgm)

(colapa <-pgpd(xx-uu,loc = 0, 
           scale = fitgm$mle[1],shape =fitgm$mle[2], lower.tail = FALSE))
(colaemp <-fitgm$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
(ret.STD2<-1/estim)
c(ret.STD2,ret.STD1)


##Con DGVE 

n <-summarise(SL, uno = sum(uno))
(n<- dplyr::select(n,uno) )
summary(n)
(n <- floor(unlist(unname(apply(n, 2,mean)))))
(mle<- fit.MSTD$mle)
(dist<-pgev(100,loc = mle[1],scale=mle[2],
           shape= mle[3]))

(estim <- 1 - dist^(1/n))
(ret.STD3<-1/estim)

###Primer periodo por arriba de 50 GMSL
m<-which.max(gmsl > 50)
kk<-length(gmsl)
mean(gmsl[m:kk]> 50)
###########Pase al siguiente script por favor.
