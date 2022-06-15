rm(list = ls())

require(SpatialExtremes)
require(latex2exp)
library(dplyr)
require(ismev)
require(cluster)
require(factoextra)
require(secr)

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

recl<-read.csv(file.choose())
names(recl)

## Datos en miles de pesos
Imp<-recl$total_claim_amount/1000
plotR(1:length(Imp),Imp, pch =20,main = "Datos")
###Dividiendo datos a ojo de buen cubero
abline(h = 2*1e1, col = "red", lwd = 2)

Km1<-which(Imp>2*1e1)
Km2<-which(Imp<=2*1e1)
Imp1<-Imp[Km1]
Imp2<-Imp[Km2] 
ylim<-c(min(Imp),max(Imp)) 

#########Datos CLusterisados
plotR(Km1,Imp1, pch =20, col = "blue"
      ,ylim = ylim,main = "Cluster Data")
points(Km2,Imp2, col = "red",pch = 20) 
###########

plotR(sort(Imp2),FEM(sort(Imp2),sort(Imp2)),type = "l",col = "purple",
      main = "FEM Datos Grupo Rojo",xlab = "u", 
      ylab = TeX("$e_{F}(u)$")) 

(u<- max(Imp1))
plotR(sort(Imp1),FEM(sort(Imp1),sort(Imp1)),type = "l",col = "purple",
      main = "FEM Datos Grupo Azul",xlab = "u", 
      ylab = TeX("$e_{F}(u)$")) 

plotR(sort(Imp),FEM(sort(Imp),sort(Imp)),type = "l",col = "purple",
      main = "FEM Todos Los Datos",xlab = "u", 
      ylab = TeX("$e_{F}(u)$")) 

########Cocientes FEM 

(n<-length(Imp2))
m<-120
plotR(sort(Imp2)[m:n],FEM(sort(Imp2),sort(Imp2))[m:n]/sort(Imp2)[m:n],type = "l",col = "purple",
      main = "Cociente FEM Datos Grupo Rojo",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$"))
####Número de observaciones por arriba de 8000 
sum(Imp2> 8000)


(n<-length(Imp1))
m<-700
plotR(sort(Imp1)[m:n],FEM(sort(Imp1),sort(Imp1))[m:n]/sort(Imp1)[m:n],type = "l",col = "purple",
      main = "Cociente FEM Datos Grupo Azul",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$")) 

n<-length(Imp)
m<-900
plotR(sort(Imp)[m:n],FEM(sort(Imp),sort(Imp))[m:n]/sort(Imp)[m:n],type = "l",col = "purple",
      main = "Cociente FEM Todos los Datos",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$")) 
#######Weibull ó Gumbel 
omegaG <- max(Imp2)
n <- length(Imp2)
Imp22 <- (omegaG - Imp2[Imp2 != omegaG])^(-1)
(n<-length(Imp22))
m<-120
summary(sort(Imp22))
plotR(sort(Imp22),FEM(sort(Imp22),sort(Imp22)),type = "l",col = "purple",
      main = "FEM Datos Grupo Rojo Transformados",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$")) 
plotR(sort(Imp22)[m:n],FEM(sort(Imp22),sort(Imp22))[m:n]/sort(Imp22)[m:n],type = "l",col = "purple",
      main = "Cociente FEM Datos Grupo Rojo Transformados",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$")) 

omegaG <- max(Imp1)
n <- length(Imp1)
Imp12 <- (omegaG - Imp1[Imp1 != omegaG])^(-1)
(n<-length(Imp12))
m<-650
summary(sort(Imp12))
plotR(sort(Imp12),FEM(sort(Imp12),sort(Imp12)),type = "l",col = "purple",
      main = "FEM Datos Grupo Azul",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$")) 
plotR(sort(Imp12)[m:n],FEM(sort(Imp12),sort(Imp12))[m:n]/sort(Imp12)[m:n],type = "l",col = "purple",
      main = "Cociente FEM Datos Grupo Azul",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$"))

###Todos los datos 
omegaG <- max(Imp)
n <- length(Imp)
Imp33 <- (omegaG - Imp[Imp != omegaG])^(-1)
(n<-length(Imp33))
m<-910
summary(sort(Imp33))
plotR(sort(Imp33),FEM(sort(Imp33),sort(Imp33)),type = "l",col = "purple",
      main = "FEM Todos los Datos",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$")) 
plotR(sort(Imp33)[m:n],FEM(sort(Imp33),sort(Imp33))[m:n]/sort(Imp33)[m:n],type = "l",col = "purple",
      main = "Cociente FEM Todos los Datos Transformados",xlab = "u", 
      ylab = TeX("$e_{F}(u)/u$"))

recl<-read.csv(file.choose())
recl <-as_tibble(recl) 
names(recl)
recl<-select(recl,c("incident_date","total_claim_amount"))
recl<-mutate(recl, incident_date =as.factor(incident_date),total_claim_amount = total_claim_amount/1000)
recl<-group_by(recl, incident_date) 
recl

tc<-summarise(recl, tc = max(total_claim_amount))$tc 

(fit1<-gev.fit(tc)) 
fit1$mle
gev.diag(fit1)

(fit2<-gum.fit(tc)) 
gum.diag(fit2)

##LogLik 

(Lik1<- exp(-fit1$nllh))
(Lik2<- exp(-fit2$nllh))

(Ratio<- 2*log(Lik1/Lik2))
###Tres parámetros en el completo, dos en el reducido
pchisq(Ratio,df = 3 - 2,lower.tail = FALSE) 
###Podemos rechazar 

###############MODIFICACIONES######################33
###Gráfica pa' seleccionar el umbral 
##No podemos elegir el umbral con el criterio porque ya no podriamos estimar
## el valor que queremos por eso de que estimos la cola en x + u
##Cola empirica en 71
mean(Imp1>71)

##Por ende, se busco en el rango de 60 a 71, i.e valores que nos permitirian 
##estimar lo que buscamos 
gpd.fitrange(Imp1,60,71) 
### Se elige 65 ? Porque se ve chidori 
u<-65
fit3<-gpd.fit(Imp1,65)
fit3$mle
##Cola empirica de G en 65
mean(Imp1>65)

###Las de diagnóstico se miran muuuy bien
gpd.diag(fit3)
##estimción del parámetro p con las Bernoulli asociadas a los datos observados
(p.hat<-mean(Imp <= 2*1e1)) 
##############FIN MODS ############################

###estim cola en 72
(x<- 72 - u)  
(colaEm <- mean(Imp1 > u)) 
(colaPare<- pgpd(x, loc = 0, scale = fit3$mle[1], shape = fit3$mle[2], lower.tail = FALSE))

###Aproximación 
(Appcola72 <- (1 - p.hat)*colaEm*colaPare)


##Por independencia
##Otra Modificación antes calculaba otra cosa
(p<-pbinom(249,500,Appcola72,lower.tail = FALSE))
