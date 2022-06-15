###Tarea extremos
###Cool plots 
rm(list =ls())
library(dplyr)
require(ismev); require(SpatialExtremes); require(latex2exp)
require(fitdistrplus);require(texmex)
require(reshape2)


##########  Cool Plots  ###########
plotR<-function(x,y,GRID = 1,...){
plot(x,y,...)
if(GRID == 1){
grid()
}
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2) 
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

##MLE of quantiles caso xi> -0.5 
VarAux<-function(alp,u,xi,a,rate){
     ifelse(xi != 0,
      u + a*(((1 - alp)/rate)^(-xi) - 1)/xi, 
      u - a*log((1-alp)/(rate)))
     }

Var<-Vectorize(VarAux,vectorize.args ="alp")

## Estimación puntual de cuantil
cuantil.ep <- function(p, muestra){
  aviso <- ""
  n <- length(muestra)
  if (p < 1/(n + 1)){
    cuantil <- min(muestra)
    aviso <- "cuantil fuera de rango muestral"
  }
  if (p > n/(n + 1)){
    cuantil <- max(muestra)
    aviso <- "cuantil fuera de rango muestral"
  }
  if ((1/(n + 1) <= p)&(p <= n/(n + 1))){
    if(((n + 1)*p)%in%(1:n)){
      j <- as.integer((n + 1)*p)
      cuantil <- sort(muestra)[j]
    }
    else {
      j <- floor((n + 1)*p)
      x.ord <- sort(muestra)[j:(j + 1)]
      cuantil <- (j + 1 - (n + 1)*p)*x.ord[1] + ((n + 1)*p - j)*x.ord[2]
    }
  }
  if (aviso != "") warning(aviso)
  return(cuantil)
}
cuantil <- function(p, x.obs) sapply(p, cuantil.ep, muestra = x.obs)

###################################

##cargar el archivo max.tex
hot <- read.table(file.choose(), header = TRUE, 
                  sep = "",dec = ".",na.strings = "M") 
set.seed(123)
for(i in 1:13){
   hot[,i] <- jitter(hot[,i])
}
hot<-hot[,-c(1,14)]
hot <- as_tibble(hot)
head(hot,3)
tail(hot,3) 
x <- melt(hot)
plt <- ggplot(data = x, aes(x = variable, y = value))
plt + geom_boxplot() + theme_minimal() + labs(x = "BoxPlot Temperatura Máximas"
      , y = "x")

##Temperaturas mínimas y máximas para temperaturas máximas
round(apply(hot,2,max,na.rm=1),3)
round(apply(hot,2,min,na.rm=1),3)
##

par(mfrow = c(3,4))
for(j in names(hot)){
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon <- sort(mon)
fem<- FEM(mon, mon)
plotR(mon, fem,col = "purple",
      main = paste("Fep MaxMax mes ", j),
      type = "l",lwd = 2) 
}

par(mfrow = c(3,4))
for(j in names(hot)){
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon <- sort(mon)
fem<- FEM(mon, mon)
plotR(mon, fem/mon,col = "purple",
      main = paste("Fep MaxMax Coc mes ", j),
      type = "l",lwd = 2) 
}

## 
u<-55
sapply(names(hot), function(j){
mon<<- dplyr::select(hot,j)
mon<<- mon[!is.na(mon)]
sum(mon > u )
}
)
x1<-c("Jan","Feb","Dec")
## 
u<-65
sapply(names(hot), function(j){
mon<<- dplyr::select(hot,j)
mon<<- mon[!is.na(mon)]
sum(mon > u )
}
)
x2<- c("Nov","Mar")
## 
u<-75
sapply(names(hot), function(j){
mon<<- dplyr::select(hot,j)
mon<<- mon[!is.na(mon)]
sum(mon > u )
}
)
x3<- c("Apr","Oct")
## 
u<-90
sapply(names(hot), function(j){
mon<<- dplyr::select(hot,j)
mon<<- mon[!is.na(mon)]
sum(mon > u )
}
)
x4<- c("May,Jun,Aug,Sep")


par(mfrow = c(3,4))
for(j in names(hot)){
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon <- sort(mon)
MM<- max(mon)
fem<- FEM(mon, mon)
if(j %in% x1){
      plotR(mon, fem/mon,col = "purple",
      main = paste("Fep MaxMax Coc. mes ", j),
      type = "l",lwd = 2, xlim = c(55,MM -1))
} 
else if(j %in% x2){
      plotR(mon, fem/mon,col = "purple",
      main = paste("Fep MaxMax Coc. mes ", j),
      type = "l",lwd = 2, xlim = c(65,MM -1))
} 
else if(j %in% x3){
      plotR(mon, fem/mon,col = "purple",
      main = paste("Fep MaxMax Coc. mes ", j),
      type = "l",lwd = 2, xlim = c(75,MM -1))
} 
else if(j %in% x4){
      plotR(mon, fem/mon,col = "purple",
      main = paste("Fep MaxMax Coc. mes ", j),
      type = "l",lwd = 2, xlim = c(90,MM -1))
} 
else{
     plotR(mon, fem/mon,col = "purple",
     main = paste("Fep MaxMax Coc. mes ", j),
     type = "l",lwd = 2, xlim = c(95,MM -1))
} 
}

####Datos transformados 

##Fep's transformados 
par(mfrow = c(3,4))
for(j in names(hot)){
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
fem<- FEM(Mon, Mon)
plotR(Mon, fem,col = "purple",
      main = paste("Fep Trans MaxMax mes ", j),
      type = "l",lwd = 2) 
}

##Cocientes Fep's transformados 
par(mfrow = c(3,4))
for(j in names(hot)){
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
fem<- FEM(Mon, Mon)
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2) 
}


##Cocientes Feps 

par(mfrow = c(2,2))
j <- "Jan"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 

j <- "Feb"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
###2by2
j <- "Mar"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
j <- "Apr"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
##2by2 *
j <- "May"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
j <- "Jun"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
## 2by2 * 
j <- "Jul"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
j <- "Aug"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
## 2 by 2  Oct *
j <- "Sep"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
j <- "Oct"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
## 
j <- "Nov"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.93)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 
j <- "Dec"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
MM<- max(mon) 
Mon <-1/(MM - mon[mon != MM])
Mon <- sort(Mon)
u<-quantile(Mon, 0.20)
q<-quantile(Mon, 0.92)
mean(Mon > q) 
fem<- FEM(Mon,Mon)
(xlim <- c(max( c(min(Mon),u)),q)) 
plotR(Mon,fem/Mon,col = "purple",
      main = paste("Fep Trans. MaxMax Coc. Mes ", j),
      type = "l",lwd = 2, xlim = xlim) 

### Analisis DGP 

j <- "Jan"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.90))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
u.jan<-55.3 
fit.jan<-gpd.fit(mon,u.jan)
round(fit.jan$mle,3)
gpd.diag(fit.jan)

j <- "Feb"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
u.feb<-55.4
fit.feb<-gpd.fit(mon,u.feb)
round(fit.feb$mle,3)
gpd.diag(fit.feb)

j <- "Mar"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
u.mar<-71.3
fit.mar<-gpd.fit(mon,u.mar)
round(fit.mar$mle,3)
gpd.diag(fit.mar)

###*
j <- "Apr"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
#no había ajuste en el primer rango propuesto 
gpd.fitrange(mon,u,q)
mean(mon>82.5); mean(mon > 83.5)
u.apr<-84.9
fit.apr<-gpd.fit(mon,u.apr)
round(fit.apr$mle,3)
gpd.diag(fit.apr)

###1* y media
j <- "May"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
u.may<-90.8
fit.may<-gpd.fit(mon,u.may)
round(fit.may$mle,3)
gpd.diag(fit.may)

###*
j <- "Jun"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
u.jun<-94.8
fit.jun<-gpd.fit(mon,u.jun)
round(fit.jun$mle,3)
gpd.diag(fit.jun)

###*
j <- "Jul"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
##Zoom para elegir 
gpd.fitrange(mon,96.6,q)
u.jul<-96.89
fit.jul<-gpd.fit(mon,u.jul)
round(fit.jul$mle,3)
gpd.diag(fit.jul)

###
j <- "Aug"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
##Zoom para elegir 
gpd.fitrange(mon,u,94.2)
u.aug<-94
fit.aug<-gpd.fit(mon,u.aug)
round(fit.aug$mle,3)
gpd.diag(fit.aug)


## *Por pp
j <- "Sep"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
##Zoom para elegir 
gpd.fitrange(mon,90,90.4)
u.sep<-90.25
fit.sep<-gpd.fit(mon,u.sep)
round(fit.sep$mle,3)
gpd.diag(fit.sep)


## *Por pp
j <- "Oct"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
##Zoom pa' elegir
gpd.fitrange(mon,82,82.5)
u.oct<-82.21
fit.oct<-gpd.fit(mon,u.oct)
round(fit.oct$mle,3)
gpd.diag(fit.oct)

## 
j <- "Nov"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
##Zoom pa' elegir
gpd.fitrange(mon,70.1,70.5)
u.nov<-70.3
fit.nov<-gpd.fit(mon,u.nov)
round(fit.nov$mle,3)
gpd.diag(fit.nov)

## *
j <- "Dec"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
(u<-quantile(mon,0.70))
(q<-quantile(mon,0.80))
mean( mon>u);mean(mon> q)
gpd.fitrange(mon,u,q)
u.dec<-59.8
fit.dec<-gpd.fit(mon,u.dec)
round(fit.dec$mle,3)
gpd.diag(fit.dec)

########## Ajustes DGVE. 

j <- "Jan"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
mle.jan <- fit$mle 

j <- "Feb"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
mle.feb <- fit$mle 

j <- "Mar"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
mle.mar <- fit$mle 

#*
j <- "Apr"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
gum.diag(gum.fit(mon))
mle.apr <- fit$mle 

#*
j <- "May"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
gum.diag(gum.fit(mon))
(mle.may <- fit$mle )

j <- "Jun"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
mle.jun <- fit$mle 

#*
j <- "Jul"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
gum.diag(gum.fit(mon))
(mle.jul <- fit$mle )

j <- "Aug"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
(mle.aug <- fit$mle)

#*
j <- "Sep"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
gum.diag(gum.fit(mon))
(mle.sep <- fit$mle)

j <- "Oct"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
gum.diag(gum.fit(mon))
(mle.oct <- fit$mle)

j <- "Nov"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
(mle.nov <- fit$mle)

j <- "Dec"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
fit<- gev.fit(mon)
gev.diag(fit)
(mle.dec <- fit$mle)

###Extremos derechos 
par(mfrow = c(2,2))
j <- "Jan"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.jan<-mean(w_F)
abline(h =w.jan,col= "red")

j <- "Feb"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.feb<-mean(w_F)
abline(h =w.feb,col= "red")

####
j <- "Mar"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.mar<-mean(w_F)
abline(h =w.mar,col= "red")

j <- "Apr"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.apr<-mean(w_F)
abline(h =w.apr,col= "red")

###May no justificado
j <- "May"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.may<-mean(w_F)
abline(h =w.may,col= "red")

j <- "Jun"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.jun<-mean(w_F)
abline(h =w.jun,col= "red")

###Jul

j <- "Jul"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.jul<-mean(w_F)
abline(h =w.jul,col= "red")

j <- "Aug"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.aug<-mean(w_F)
abline(h =w.aug,col= "red")

####

j <- "Sep"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.sep<-mean(w_F)
abline(h =w.sep,col= "red")

j <- "Oct"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.oct<-mean(w_F)
abline(h =w.oct,col= "red")

####

j <- "Nov"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.nov<-mean(w_F)
abline(h =w.nov,col= "red")

j <- "Dec"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)]
mon<-sort(mon)
n<-length(mon)
w_F <- ED_estimator(mon,1)
for (k in 2:(floor(n^0.8))){
  w_F <- c(w_F,ED_estimator(mon,k))
}
plotR(1:length(w_F),w_F, type="l", col="purple",
     main = paste("Extremo derecho ",j),  lwd = 2)
w.dec<-mean(w_F)
abline(h =w.dec,col= "red")

####etstremos estimados 
round(c(w.jan,w.feb,w.mar,w.apr,w.may,w.jun,
  w.jul,w.aug,w.sep,w.oct,w.nov,w.dec),3)

####
htn <- names(hot)[names(hot) != "May"]
par(mfrow = c(3,4))
for (j in htn){
	mon<- dplyr::select(hot,j)
	mon <- mon[!is.na(mon)]
	mon<-sort(mon)
	n<-length(mon)
	p_negative <- test_zero_vs_negative(mon,10)
      for (k in 11:(floor(n^0.8))){
  		p_negative <- c(p_negative, test_zero_vs_negative(mon,k))
	}
	plotR(10+1:length(p_negative),p_negative, type = "l"
      ,col = "purple",lwd  =2,
      main = paste("Prueba de hip, dom de atrac. ",j))
	abline(h = 0.05, col = "red")
}
###

######Excesos
j <- "Jan"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.jan - mon)
summary(mon)
uu <-1/(w.jan - u.jan)
xx <-1/(w.jan - MM)
fit1<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit1$mle[1],shape =fit1$mle[2], lower.tail = FALSE))
(colaemp <-fit1$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

j <- "Feb"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.feb - mon)
summary(mon)
uu <-1/(w.feb - u.feb)
xx <-1/(w.feb - MM)
fit2<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit2$mle[1],shape =fit2$mle[2], lower.tail = FALSE))
(colaemp <-fit2$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

j <- "Mar"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.mar - mon)
summary(mon)
uu <-1/(w.mar - u.mar)
xx <-1/(w.mar - MM)
fit3<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit3$mle[1],shape =fit3$mle[2], lower.tail = FALSE))
(colaemp <-fit3$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

j <- "Apr"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.apr - mon)
summary(mon)
uu <-1/(w.apr - u.apr)
xx <-1/(w.apr - MM)
fit4<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit4$mle[1],shape =fit4$mle[2], lower.tail = FALSE))
(colaemp <-fit4$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

j <- "May"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.may- mon)
summary(mon)
uu <-1/(w.may - u.may)
xx <-1/(w.may - MM)
fit5<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit5$mle[1],shape =fit5$mle[2], lower.tail = FALSE))
(colaemp <-fit5$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
1/estim

j <- "Jun"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.jun- mon)
summary(mon)
uu <-1/(w.jun - u.jun)
xx <-1/(w.jun - MM)
fit6<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit6$mle[1],shape =fit6$mle[2], lower.tail = FALSE))
(colaemp <-fit6$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

j <- "Jul"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.jul- mon)
summary(mon)
uu <-1/(w.jul - u.jul)
xx <-1/(w.jul - MM)
fit7<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit7$mle[1],shape =fit7$mle[2], lower.tail = FALSE))
(colaemp <-fit7$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

j <- "Aug"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.aug- mon)
summary(mon)
uu <-1/(w.aug - u.aug)
xx <-1/(w.aug - MM)
fit8<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit8$mle[1],shape =fit8$mle[2], lower.tail = FALSE))
(colaemp <-fit8$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

j <- "Sep"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.sep- mon)
summary(mon)
uu <-1/(w.sep - u.sep)
xx <-1/(w.sep - MM)
fit9<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit9$mle[1],shape =fit9$mle[2], lower.tail = FALSE))
(colaemp <-fit9$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)


j <- "Oct"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.oct- mon)
summary(mon)
uu <-1/(w.oct- u.oct)
xx <-1/(w.oct - MM)
fit10<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit10$mle[1],shape =fit10$mle[2], lower.tail = FALSE))
(colaemp <-fit10$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

j <- "Nov"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.nov- mon)
summary(mon)
uu <-1/(w.nov- u.nov)
xx <-1/(w.nov - MM)
fit11<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit11$mle[1],shape =fit11$mle[2], lower.tail = FALSE))
(colaemp <-fit11$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

j <- "Dec"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
Mon <-1/(w.dec- mon)
summary(mon)
uu <-1/(w.dec- u.dec)
xx <-1/(w.dec - MM)
fit12<-gpd.fit(Mon,uu)
(colapa <-SpatialExtremes::pgpd(xx -uu,loc = 0, 
           scale = fit12$mle[1],shape =fit12$mle[2], lower.tail = FALSE))
(colaemp <-fit12$rate)
(estim<- colapa*colaemp)
###Periodo de retorno
round(1/estim,3)

#####Sin excesos
j <- "Jan"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.jan[1],scale=mle.jan[2],
                      shape= mle.jan[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "Feb"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.feb[1],scale=mle.feb[2],
                      shape= mle.feb[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "Mar"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.mar[1],scale=mle.mar[2],
                      shape= mle.mar[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "Apr"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.apr[1],scale=mle.apr[2],
                      shape= mle.apr[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "May"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.may[1],scale=mle.may[2],
                      shape= mle.may[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "Jun"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.jun[1],scale=mle.jun[2],
                      shape= mle.jun[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "Jul"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.jul[1],scale=mle.jul[2],
                      shape= mle.jul[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "Aug"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.aug[1],scale=mle.aug[2],
                      shape= mle.aug[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "Sep"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.sep[1],scale=mle.sep[2],
                      shape= mle.sep[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "Oct"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.oct[1],scale=mle.oct[2],
                      shape= mle.oct[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

j <- "Nov"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
(ee<-SpatialExtremes::pgev(MM,loc = mle.nov[1],scale=mle.nov[2],
                      shape= mle.nov[3],lower.tail =  FALSE))
##Periodo de retorno GEV
round(1/ee,3)

j <- "Dec"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
MM<-max(mon)
ee<-SpatialExtremes::pgev(MM,loc = mle.dec[1],scale=mle.dec[2],
                      shape= mle.dec[3],lower.tail =  FALSE)
##Periodo de retorno GEV
round(1/ee,3)

##################VAR'S 
#####################
##Nov
ff<- fit.nov
(u<- ff$threshold)
(rate<- ff$rate)
(a<- ff$mle[1])
(xi<- ff$mle[2])
alp <- c(0.95,0.99)

####Var Para Nov:
if(xi >-0.5){      
    Var(alp,u,xi,a,rate)
}
###Data trans otra forma de estimar, que esta aprueba
ff<- fit11
u<- ff$threshold
rate<- ff$rate
a<- ff$mle[1]
xi<- ff$mle[2]
w.nov - 1/Var(alp,u,xi,a,rate)

j <- "Nov"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
##Método
cuantil(alp, mon)

##Dec
ff<- fit.dec
(u<- ff$threshold)
(rate<- ff$rate)
(a<- ff$mle[1])
(xi<- ff$mle[2])
alp <- c(0.95,0.99)

####Var Para Dec:
if(xi >-0.5){      
    Var(alp,u,xi,a,rate)
}
###Data trans otra forma de estimar, que esta aprueba
ff<- fit12
u<- ff$threshold
rate<- ff$rate
a<- ff$mle[1]
xi<- ff$mle[2]
w.dec - 1/Var(alp,u,xi,a,rate)

j <- "Dec"
mon<- dplyr::select(hot,j)
mon <- mon[!is.na(mon)] 
##Método
cuantil(alp, mon)

######Estimaciones 
######TLC aplica en cada caso porque las temps no pesentan evidencia en contra
######De extremo derecho finito. 
(a<-apply(hot,2,mean,na.rm = 1))


u<-c(
     fit.jan$rate,fit.feb$rate,fit.jan$rate,fit.mar$rate
     ,fit.apr$rate,fit.may$rate,fit.jun$rate,fit.jul$rate,fit.aug$rate
     ,fit.sep$rate,fit.oct$rate,fit.nov$rate,fit.dec$rate
     )

(n<-sapply(names(hot),function(j)sum(!is.na(dplyr::select(hot,j)))))
s<-length(u)
sapply(1:s,function(j)(1 - u[j])^(n[j]))

