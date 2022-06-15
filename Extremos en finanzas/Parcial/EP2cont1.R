######################### Análisis diferencias GMSL 
 j<- "gmsl"
gmsl<- unname(unlist(dplyr::select(SeaLev,all_of(j))))
xxxx<-gmsl[1]

(gmsl<-c(gmsl[1],diff(gmsl)))


par(mfrow = c(1,1))
plotR(gmsl, type = "p", pch = 20, col = "red", 
       main = paste("Gráfico de dispersión: Diff. GMSL"),
       ylab = "y")

summary(gmsl)
####### Inciso d) Determinemos si existe un dominio de atracción máximal 
#######
###FEP's 

par(mfrow = c(1,2))
gmsl <-sort(gmsl)
ff1<-FEM(gmsl, gmsl)
plotR(gmsl,ff1, type = "l", col = "purple", lwd = 2,
       main = "FEP Diff GMSL")
abline(v = 5, col = "red")
sum(gmsl > 5)
u<-quantile(gmsl,0.65)
q<-quantile(gmsl,0.985)
xlim = c(u,q)
plotR(gmsl,ff1, type = "l", col = "purple", lwd = 2,
       main = "FEP Diff GMSL",xlim = xlim,ylim = c(1.3,2))
mean(gmsl > 5)

##Coc Fep 
gmsl <-sort(gmsl)
ff1<-FEM(gmsl, gmsl)
plotR(gmsl,ff1/gmsl, type = "l", col = "purple", lwd = 2,
      main = "Cociente FEP empír.Diff GMSL")
u<-quantile(gmsl,0.65)
q<-quantile(gmsl,0.985)
mean(gmsl > u);mean(gmsl > q)
xlim<-c(u,q) 
plotR(gmsl,ff1/gmsl, type = "l", col = "purple", lwd = 2
      ,xlim=xlim,ylim = c(0,1.5),  main = "Cociente FEP empír. Diff GMSL Zoom")
)

######################### Datos transformados ######################: 

######################################################Ajustes DGP :DD

(u<-quantile(gmsl,0.95))
(q<-quantile(gmsl,0.94))
mean(gmsl > u);mean(gmsl > q)
gpd.fitrange(gmsl,q,u) 
###Zoom elección umbraaal 
gpd.fitrange(gmsl,q,4.15) 
u.gmsl<-4.10
fit.gmsl<-gpd.fit(gmsl,u.gmsl)
gpd.diag(fit.gmsl)
round(fit.gmsl$mle,3) 
#############Prueba G_{n,k}(0)
mon<-sort(gmsl)
n<-length(mon)
p_negative <- test_zero_vs_negative(mon,1)
      for (k in 2:(floor(n^0.7))){
  		p_negative <- c(p_negative, test_zero_vs_negative(mon,k))
	}

p_positive <- test_zero_vs_positive(mon,1)
      for (k in 2:(floor(n^0.7))){
  		p_positive <- c(p_positive, test_zero_vs_positive(mon,k))
	}
par(mfrow = c(1,2))
plotR(p_negative,type = "l",lwd = 2, col ="purple",
      main = TeX("STD $\\xi < 0$ vs $\\xi = 0$"))
plotR(p_positive,type = "l",lwd = 2, col ="purple",
      main = TeX("STD $\\xi > 0$ vs $\\xi = 0$"))
abline(h = 0.05, col = "red")

##########################################################
(x<-xxxx)
aaa<- 0.5
(a<-fit.gmsl$mle[1])
(xi <-fit.gmsl$mle[2])
(fu<-fit.gmsl$rate)
sl<-summarise(SL, ss = sum(uno))
sl <- data.frame(sl)
mean(sl$ss[-length(sl$ss)])
sum(sl$ss)
(n<-1048 + 19 + 3*37)
#Estim
x+u.gmsl +a*(((1-aaa)/(n*fu))^(-xi) - 1)/xi 

############################################## Proba solicitada
######################################################Bloques de tamaño 8


mgmslm <- sapply(1:131,function(j)max(gmsl[(8*(j-1) +1):(8*j)])) 
mgmslm
(xx<-max(gmsl))
fit.g2<-gev.fit(mgmslm)
gev.diag(fit.g2)
round(fit.g2$mle,3) 
(mle<- fit.g2$mle)

###Probabilidad.
pgev(xx,loc = mle[1],scale = mle[2],shape =mle[3],lower.tail =FALSE)

### \PP[diff(X) > xx]\overline{F}(xx) con a como abajo
### Por método de excesos:
### \overline{F}(xx) \approx \overline{F}_{n}(u.gmsl)\overline{Patito}(xx - u.gmsl)

(xx<-max(gmsl))
(colapa <-pgpd(xx-u.gmsl,loc = 0, 
           scale = fit.gmsl$mle[1],shape =fit.gmsl$mle[2], lower.tail = FALSE))
(colaemp <-fit.gmsl$rate)
colapa*colaemp
(estim<- 1-(1-colapa*colaemp)^8)


#############################Fep 
#####  E[diff(X) | diff(X) > 2.96]= E[diff(X) -2.96| diff(X) > 2.96] + 2.96
#####                             = e_{F}(2.96) + 2.96
FEM(2.96,gmsl) + 2.96
(alp<-1/fit.gmsl$mle[2])
2.96*(alp/(alp - 1))
 
mean(gmsl>2.96)
quan

############# continued...
