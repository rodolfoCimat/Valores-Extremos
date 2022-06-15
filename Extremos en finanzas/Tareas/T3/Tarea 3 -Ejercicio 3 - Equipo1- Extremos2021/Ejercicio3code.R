rm(list = ls())

require(SpatialExtremes)
require(latex2exp)
library(dplyr)
require(ismev)
require(cluster)
require(factoextra)
require(secr)
require(fitdistrplus)

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
###################################
data("danishuni")
data <- danishuni
names(data)
Imp<-data$Loss
plotR(1:length(Imp),Imp, pch =20,main = "Datos")

plotR(sort(Imp),FEM(sort(Imp),sort(Imp)),type = "l",col = "purple",
      main = "FEM Todos Los Datos",xlab = "u", 
      ylab = TeX("$e_{F}(u)$")) 

(n<-length(Imp))
m<-900
plotR(sort(Imp)[m:n],FEM(sort(Imp),sort(Imp))[m:n]/sort(Imp)[m:n],type = "l",col = "purple",
      main = "Cociente FEM Todos los Datos",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$")) 

(n<-length(Imp)-300)
m<-900
plotR(sort(Imp)[m:n],FEM(sort(Imp),sort(Imp))[m:n]/sort(Imp)[m:n],type = "l",col = "purple",
      main = "Cociente FEM Todos los Datos",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$")) 


data1<-as_tibble(cbind(data,rep(1, nrow(data))))
(data1<-group_by(data1, Date))
names(data1)[3]<-"uno"
###reclamos máximos por día 
(tc<-summarise(data1, tc = max(Loss)))
head(tc,3)
tail(tc,3)
tc<-tc$tc

###reclamos acumulados por día 
(sc<-summarise(data1, sc = sum(uno)))
###Como puede ver hay días sin reclamaciones, 

###Se rrelenaran con ceros dichos dias para representar que en los mismos
###No hubo reclamaciones
Date1<-seq(from =as.Date("1980-01-03"), to =as.Date("1990-12-31"),by = 1  )

sc<-sapply(Date1, function(x)
     as.numeric(ifelse(x %in% sc$Date,sc[sc$Date == x,2], 0))
     )
sc<-as_tibble(data.frame(date = Date1, n.claims = sc) )
head(sc,3)
tail(sc,3)
sc<-sc$n.claims

####Ajuste generalizada de extremos a máximos por dias
(fit1<-gev.fit(tc)) 
gev.diag(fit1)
fit1$mle

##Empírica en x = 9.5
mean(Imp> 9.5)
sum(Imp>9.5)

gpd.fitrange(Imp,8,9.5) 
### Se elige 8.25 
u<-8.25
###
mean(Imp> 8.25)
sum(Imp>8.25)

###Ajuste pareto
fit3<-gpd.fit(Imp,12)
(l.hat<-mean(sc))

###Diagnóstico de ajuste 
gpd.diag(fit3) 

###estim cola en 10
(x<- 10 - u)  
(colaEm <- mean(Imp > u)) 
(colaPare<- pgpd(x, loc = 0, scale = fit3$mle[1], shape = fit3$mle[2], lower.tail = FALSE))

###Aproximación
##Número de observaciones Poisson(lambda) independientes 
(kplusm<-length(sc))
##Número de observaciones de los montos de reclamación 
(m<-nrow(data1))
(k<-kplusm - m)

t <- seq(1e-6,210, length.out = 1000)
Nt1<-t*l.hat*colaEm*colaPare
Nt2<- 1-exp(-l.hat*t*colaEm*colaPare)
Nt1[Nt1 > 1]<- 1
plotR(t,Nt1,col ="red", type = "l",lwd = 2,
      main = TeX("Approx. Cola de $M_{N(t)}$ Evaluada en x = 10"),
      xlab = "t")
lines(t,Nt2, col = "blue", lwd = 2)
