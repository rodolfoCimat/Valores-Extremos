#######Diferencias STD 

j<- "STD"
(STD<- unname(unlist(dplyr::select(SeaLev,all_of(j)))))
(STD<-c(STD[1],diff(STD)))
(std1<- STD) 
par(mfrow = c(1,1))
plotR(STD, type = "p", pch = 20, col = "red", 
       main = paste("Gráfico de dispersión: Diff STD"),ylab ="y")

summary(STD)
####### FEP 
STD <-sort(STD)
ff1<-FEM(STD,STD)
par(mfrow = c(1,2))
plotR(STD,ff1, type = "l", col = "purple", lwd = 2,
      main = "FEP Diff STD")###0.010
abline(v = 93)##0.084
u<-quantile(STD,0.65)
q<-quantile(STD,0.985)
xlim = c(u,q)
plotR(STD,ff1, type = "l", col = "purple", lwd = 2,
       main = "FEP Diff STD",xlim = xlim,ylim = c(1,7))
mean(gmsl > 5)

###### FEP cocientada No evidencia vs Frechet Domain 
ff2<-FEM(STD, STD)
plotR(STD,ff2/STD, type = "l", col = "purple", lwd = 2,
      main = "Cociente FEP empír. Diff STD")
(u<-quantile(STD,0.65))
(q<-quantile(STD,0.985))
mean(STD > u);mean(STD > q)
xlim<-c(u,q) 
plotR(STD,ff2/STD,type = "l", col = "purple", lwd = 2
      ,xlim=xlim,main = "Cociente FEP empír Diff STD Zoom",
       ylim = c(1,2))
###################################################### Pruebas G_{n,k}(0)

mon<-sort(STD)
n<-length(mon)
p_negative <- test_zero_vs_negative(mon,1)
      for (k in 2:(floor(n^0.8))){
  		p_negative <- c(p_negative, test_zero_vs_negative(mon,k))
	}

p_positive <- test_zero_vs_positive(mon,1)
      for (k in 2:(floor(n^0.8))){
  		p_positive <- c(p_positive, test_zero_vs_positive(mon,k))
	}
summary(p_negative)
plotR(p_negative,type = "l",lwd = 2, col ="purple")
abline(h = 0.05, col = "red")
plotR(p_positive,type = "l",lwd = 2, col ="purple")

######################################################Ajustes DGP :DD

(u<-quantile(STD,0.95))
(q<-quantile(STD,0.94))
mean(STD > u);mean(STD > q)
gpd.fitrange(STD,q,u) 
###Zoom elección umbraaal  
gpd.fitrange(STD,2.80,2.86)
u.STD<-2.83
fit.STD<-gpd.fit(STD,u.STD)
gpd.diag(fit.STD)

##############################################Estimación mediante la mediana

(x<-max(STD))
aaa<- 0.5
a<-fit.STD$mle[1]; xi <-fit.STD$mle[2]
fu<-fit.STD$rate
sl<-summarise(SL, ss = sum(uno))
sl <- data.frame(sl)
mean(sl$ss[-length(sl$ss)])
sum(sl$ss)
n<-1048 + 19 + 3*37
#Estim
x+u.STD +a*(((1-aaa)/(n*fu))^(-xi) - 1)/xi 

############################################### DGVE
mgmslm <- sapply(1:131,function(j)max(std1[(8*(j-1) +1):(8*j)])) 
mgmslm
(xx<-c(max(STD),max(STD[STD!= max(STD)])))
fit.g2<-gev.fit(mgmslm)
gev.diag(fit.g2)
round(fit.g2$mle,3) 
(mle<- fit.g2$mle)
pgev(xx,loc = mle[1],scale = mle[2],shape =mle[3],lower.tail =FALSE)

############################################## Proba solicitada
### \PP[diff(X) > xx]\overline{F}(xx) con a como abajo
### Por método de excesos:
### \overline{F}(xx) \approx \overline{F}_{n}(u.gmsl)\overline{Patito}(xx - u.gmsl)

(colapa <-pgpd(xx-u.STD,loc = 0, 
           scale = fit.STD$mle[1],shape =fit.STD$mle[2], lower.tail = FALSE))
(colaemp <-fit.STD$rate)
(estim<- colapa*colaemp)
(estim<- 1-(1-colapa*colaemp)^8)
###Periodo de retorno
(ret.STD<-1/estim)

#############################Fep 
#####  E[diff(X) | diff(X) > 2.96]= E[diff(X) -2.96| diff(X) > 2.96] + 2.96
#####                             = e_{F}(2.96) + 2.96
FEM(2.96,STD) + 2.96
(alp<-1/fit.STD$mle[2])
2.96*(alp/(alp - 1))
