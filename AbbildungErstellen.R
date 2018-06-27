library(ggplot2) #zeichnen
library(MASS)
library(grid)#Lokation von ggplots2
library(fGarch)#Skewed-t-Verteilung und GARCH
library(reshape2)
library(extRemes) #Extremwerttheorie
library(xts)      #eXtensible Time Series
library(tseries)  #Zeitreihe
library(MASS)
library(rugarch)  #univariate GARCH
library(skewt)    #Skewed-t-Verteilung
library(forecast) #Forecasting time series 
library(evir)   #auch Extremwerttheorie

#Abbildung 1: Studentsche t-Verteilung und NV-Verteilung
tfitt=fitdistr(-dax_log$logreturn,"t")
tfitt
m=tfitt$estimate[1]
s=tfitt$estimate[2]
df=tfitt$estimate[3]
pi=3.1415926
ft=function(x,df,m,s){(gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5))}

t_nv=ggplot(data=dax_log,aes(-logreturn))+ labs(x = "Log-Verlust",y="Dichte")+
geom_histogram(aes(y=..density..),fill="blue",bins=100)+
#geom_line(aes(y = ..density..), stat = 'density',lwd=1.2)+
  stat_function(fun = dnorm, n = 1000, args = list(mean = mean(-dax_log$logreturn), sd = sd(-dax_log$logreturn)),lwd=1.1) +
  stat_function(fun = function(x) (gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5)) ,lwd=1.1,col="red") +
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),)
t_nv


t_nv1=ggplot(data=dax_log,aes(-logreturn))+ labs(x = "Log-Verlust",y="Dichte")+coord_cartesian(xlim=c(3, 8),ylim=c(0, 0.05))+
  geom_histogram(aes(y=..density..),fill="blue",bins=200)+
  #geom_line(aes(y = ..density..), stat = 'density',lwd=1.1,col="green")+
  stat_function(fun = dnorm, n = 1000, args = list(mean = mean(-dax_log$logreturn), sd = sd(-dax_log$logreturn)),lwd=1.1) +
  stat_function(fun = function(x) (gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5)) ,lwd=1.1,col="red") +
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),)
t_nv1

##Location
vp1= viewport(width=0.5,height=1,x = 0, y = 1, just = c("left", "top"))
vp2= viewport(width=0.5,height=1,x = 0.75, y = 1, just = c("center", "top"))

print(t_nv,vp=vp1)
print(t_nv1,vp=vp2)


#Abbildung 2: Skewed t Verteilung

sktfit=sstdFit(-dax_log$logreturn)
sktfit$estimate

skt=ggplot(data=dax_log,aes(-logreturn))+ labs(x = "Log-Verlust",y="Dichte") +#coord_cartesian(xlim=c(3, 8),ylim=c(0, 0.05))+
  geom_histogram(aes(y=..density..),fill="blue",bins=200)+
  #labs(title="Histogramm des Log-Verlusts")+
  #geom_line(aes(y = ..density..), stat = 'density',lwd=1.1,col="green")+
  stat_function(fun = dnorm, n = 1000, args = list(mean = mean(-dax_log$logreturn), sd = sd(-dax_log$logreturn)),lwd=1.1) +
  stat_function(fun = function(x) (gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5)) ,lwd=1.1,col="red") +
  stat_function(fun = dskt, n = 1000, args = list(df=sktfit$estimate[3],gamma=sktfit$estimate[4]),lwd=1.1,col="green") +
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),)
skt

skt1=ggplot(data=dax_log,aes(-logreturn))+ labs(x = "Log-Verlust",y="Dichte") +coord_cartesian(xlim=c(3, 8),ylim=c(0, 0.05))+
  geom_histogram(aes(y=..density..),fill="blue",bins=200)+
  #labs(title="Histogramm des Log-Verlusts")+
  #geom_line(aes(y = ..density..), stat = 'density',lwd=1.1,col="green")+
  stat_function(fun = dnorm, n = 1000, args = list(mean = mean(-dax_log$logreturn), sd = sd(-dax_log$logreturn)),lwd=1.1) +
  stat_function(fun = function(x) (gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5)) ,lwd=1.1,col="red") +
  stat_function(fun = dskt, n = 1000, args = list(df=sktfit$estimate[3],gamma=sktfit$estimate[4]),lwd=1.1,col="green") +
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),)
skt1

print(skt,vp=vp1)
print(skt1,vp=vp2)

#Abbildung 3: EVT Dichte

evtd=ggplot(data.frame(x=c(-4, 4)), aes(x)) +labs(x = "x",y="Dichte")  +
  stat_function(fun=devd,args = list(loc = 0, scale = 1, shape = 0, threshold = 0, log = F,type="GEV"),lwd=1.1)+
  stat_function(fun=devd,args = list(loc = 1, scale = 1, shape = 1, threshold = 0, log = F,type="GEV"),lwd=1.1,col="red")+
  stat_function(fun=devd,args = list(loc = -1, scale = 1, shape = -1, threshold = 0, log = F,type="GEV"),lwd=1.1,col="green")+
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),)
evtd


#evt mit Legend

x  <- seq(-5, 5, 0.01) 
Gumbel <- devd(x,loc = 0, scale = 1, shape = 0, threshold = 0, log = F,type="GEV")
Frechet <- devd(x,loc = 1, scale = 1, shape = 1, threshold = 0, log = F,type="GEV")
Weibull <- devd(x,loc = -1, scale = 1, shape = -1, threshold = 0, log = F,type="GEV")
df <- data.frame(x,Gumbel,Frechet,Weibull)

gg <- melt(df,id="x")
ggplot(gg, aes(x=x, y=value, color=variable))+
  geom_line(lwd=1.1)+labs(x = "x",y="Dichte")+
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),legend.text=element_text(size=17),legend.title=element_blank())

#Abbildung 4: GPD


x  <- seq(0.01, 5, 0.01) 
gpd0 <- devd(x,loc = 0, scale = 1, shape = 0, threshold = 0, log = F,type="GP")
gpd1 <- devd(x,loc = 0, scale = 1, shape = 0.5, threshold = 0, log = F,type="GP")
gpd2 <- devd(x,loc = 0, scale = 1, shape = -0.5, threshold = 0, log = F,type="GP")

gpd012 <- data.frame(x,gpd0,gpd1,gpd2)

ggpd <- melt(gpd012,id="x")
gpd11=ggplot(ggpd, aes(x=x, y=value, color=factor(variable,labels=c("xi=0","xi=0.5","xi=-0.5"))))+
  geom_line(lwd=1.1)+labs(x = "y",y="g(y)")+
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),legend.text=element_text(size=17),legend.title=element_blank())

x  <- seq(0, 5, 0.01) 
pgpd0 <- pevd(x,loc = 0, scale = 1, shape = 0, threshold = 0, log = F,type="GP")
pgpd1 <- pevd(x,loc = 0, scale = 1, shape = 0.5, threshold = 0, log = F,type="GP")
pgpd2 <- pevd(x,loc = 0, scale = 1, shape = -0.5, threshold = 0, log = F,type="GP")

pgpd012 <- data.frame(x,pgpd0,pgpd1,pgpd2)

pggpd <- melt(pgpd012,id="x")
gpd22=ggplot(pggpd, aes(x=x, y=value, color=factor(variable,labels=c("xi=0","xi=0.5","xi=-0.5"))))+
  geom_line(lwd=1.1)+labs(x = "y",y="G(y)")+
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),legend.text=element_text(size=17),legend.title=element_blank())
  
print(gpd11,vp=vp1)
print(gpd22,vp=vp2)

#Abbildung 5: Mean Residual Life Plot (Mean Excess)
mrlplot(dax_log_xts, main="Sample Mean Residual Life Plot",)    #u ist vielleicht in (0,4), aber nicht informative
meplot(dax_log_xts,xlim=c(0,5),ylim=c(1,1.5),type="l",main="Sample Mean Residual Life Plot")  #u ist vielleicht 3.5. Nach 3.5 ist linear. Ist (6826-100)/6826=0.9854Quantil

#Abbildung 6: Hill-Plot
par(mfrow=c(1,2))
hill(dax_log_xts)
hill(dax_log_xts,xlim=c(60,340)) 


#Abbildung 7:Plot der Zeitreihe
plot.xts(dax_xts,main="Schlusskurs des DAX Index")

#Abbildung 8: Plot der Log-Verluste
plot(dax_log_xts,main="Log-Verlust")

#Abbildung 9: Norm-QQ-plot
stats::qqnorm(dax_log_xts);qqline(dax_log_xts)

#Abbildung 10: Monatliche Maxima
monatmax=apply.monthly(dax_log_xts,max) #die monatliche groesste Verlust 
plot(monatmax,main="Monatliche Maxima")

#Abbildung 11: VaR-Prognosen der BM-Methode
plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
legend("topright",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 12: ES-Prognosen der BM-Methode
plot(dax_log_xts[1051:6826],main="ES",,ylim=c(0,10))  
lines(ES995_bmm_xts,col="red")   
lines(ES99_bmm_xts,col="blue")    
lines(ES95_bmm_xts,col="green")  
legend("topright",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 13: VaR-Prognosen der BM-Methode mit dem extremen Index
plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_bmme_xts,col="red")   
lines(VaR99_bmme_xts,col="blue")    
lines(VaR95_bmme_xts,col="green")  
legend("topright",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 14: ES-Prognosen der BM-Methode mit dem extremen Index
plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
lines(ES995_bmme_xts,col="red")   
lines(ES99_bmme_xts,col="blue")    
lines(ES95_bmme_xts,col="green")  
legend("topright",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 15: VaR-Prognosen der POT-Methode 
plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_pot_xts,col="red")   
lines(VaR99_pot_xts,col="blue")    
lines(VaR95_pot_xts,col="green") 
legend("topright",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 16: ES-Prognosen der POT-Methode 
plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10))  
lines(ES995_pot_xts,col="red")   
lines(ES99_pot_xts,col="blue")    
lines(ES95_pot_xts,col="green") 
legend("topright",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 17: Betrag der Log-Verluste
plot(abs(dax_log_xts),main="Betrag der Log-Verluste")

#Abbildung 18: Standardabweichung der Log-Verluste
plot(garchfit3@sigma.t,type="l",main="Standardabweichung der Reihe")

#Abbildung 19: Standardisierte Residuen des GARCH(1,1)-Modells
plot(zt,main="Standardisierte Residuen")

#Abbildung 20: Autokorrelationsfunktion
par(mfrow=c(2,2))
acf(dax_log_xts,main="ACF der Log-Verluste");acf(abs(dax_log_xts),main="ACF der Abs(Log-Verluste)")  #nicht i.i.d
acf(zt,main="ACF der Standardisierten Residuen");acf(abs(zt),main="ACF der Abs(Standardisierten Residuen)")     #keinen ARCH-Effekt
par(mfrow=c(1,1))

#Abbildung 21: VaR-Prognosen der POT-Methode mit dem GARCH-Filter
plot(dax_log_xts[1051:2050],main="VaR",ylim=c(0,12))  
lines(VaR995_pot_xts_garch0[1:1000],col="red")   
lines(VaR99_pot_xts_garch0[1:1000],col="blue")    
lines(VaR95_pot_xts_garch0[1:1000],col="green")   
legend("topleft",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 22: ES-Prognosen der POT-Methode mit dem GARCH-Filter
plot(dax_log_xts[1051:2050],main="ES",ylim=c(0,13))  
lines(ES995_pot_xts_garch0[1:1000],col="red")   
lines(ES99_pot_xts_garch0[1:1000],col="blue")    
lines(ES95_pot_xts_garch0[1:1000],col="green") 
legend("topleft",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 23: VaR-Prognosen der historischen Simulation
plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_hs,dax_log$date[1051:6826]),col="red")   
lines(xts(VaR99_hs,dax_log$date[1051:6826]),col="blue")    
lines(xts(VaR95_hs,dax_log$date[1051:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 24: ES-Prognosen der historischen Simulation
plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10))  
lines(xts(ES995_hs,dax_log$date[1051:6826]),col="red")   
lines(xts(ES99_hs,dax_log$date[1051:6826]),col="blue")    
lines(xts(ES95_hs,dax_log$date[1051:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 25: VaR-Prognosen der normalen Verteilung
plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_nv,dax_log$date[1051:6826]),col="red")   
lines(xts(VaR99_nv,dax_log$date[1051:6826]),col="blue")    
lines(xts(VaR95_nv,dax_log$date[1051:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 26: ES-Prognosen der normalen Verteilung 
plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10))  
lines(xts(ES995_nv,dax_log$date[1051:6826]),col="red")   
lines(xts(ES99_nv,dax_log$date[1051:6826]),col="blue")    
lines(xts(ES95_nv,dax_log$date[1051:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 27: VaR-Prognosen der Studentschen t-Verteilung
plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_t,dax_log$date[1051:6826]),col="red")   
lines(xts(VaR99_t,dax_log$date[1051:6826]),col="blue")    
lines(xts(VaR95_t,dax_log$date[1051:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 28: ES-Prognosen der Studentschen t-Verteilung 
plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
lines(xts(ES995_t,dax_log$date[1051:6826]),col="red")   
lines(xts(ES99_t,dax_log$date[1051:6826]),col="blue")    
lines(xts(ES95_t,dax_log$date[1051:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 29: VaR-Prognosen der Skewed t-Verteilung
plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(xts(-VaR995_skt1,dax_log$date[1051:6826]),col="red")   
lines(xts(-VaR99_skt1,dax_log$date[1051:6826]),col="blue")    
lines(xts(-VaR95_skt1,dax_log$date[1051:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 30: ES-Prognosen der Skewed t-Verteilung
plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
lines(xts(-ES995_skt1,dax_log$date[1051:6826]),col="red")   
lines(xts(-ES99_skt1,dax_log$date[1051:6826]),col="blue")    
lines(xts(-ES95_skt1,dax_log$date[1051:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 31: VaR-Prognosen des GARCH-Modells mit der normalen Verteilung 
plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_garchnv_xts,col="red")   
lines(VaR99_garchnv_xts,col="blue")    
lines(VaR95_garchnv_xts,col="green")   
legend("topright",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 32: ES-Prognosen des GARCH-Modells mit der normalen Verteilung 
plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
lines(ES995_garchnv_xts,col="red")   
lines(ES99_garchnv_xts,col="blue")    
lines(ES95_garchnv_xts,col="green") 
legend("topright",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 33: VaR-Prognosen des GARCH-Modells mit der Studentschen t-Verteilung
plot(dax_log_xts[1051:2050],main="VaR",ylim=c(0,12))  
lines(VaR995_garcht_xts[1:1000],col="red")   
lines(VaR99_garcht_xts[1:1000],col="blue")    
lines(VaR95_garcht_xts[1:1000],col="green")   
legend("topleft",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 34: ES-Prognosen des GARCH-Modells mit der Studentschen t-Verteilung .
plot(dax_log_xts[1051:2050],main="ES",ylim=c(0,13))  
lines(ES995_garcht_xts[1:1000],col="red")   
lines(ES99_garcht_xts[1:1000],col="blue")    
lines(ES95_garcht_xts[1:1000],col="green") 
legend("topleft",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 35: VaR-Prognosen des GARCH-Modells mit der Skewed t-Verteilung
plot(dax_log_xts[1051:2050],main="VaR",ylim=c(0,12))  
lines(VaR995_garchsktnew_xts[1:1000],col="red")   
lines(VaR99_garchsktnew_xts[1:1000],col="blue")    
lines(VaR95_garchsktnew_xts[1:1000],col="green")   
legend("topleft",inset=0.005,c("VaR0.005","VaR0.01","VaR0.05"),col=c("red","blue","green"),lty=c(1,1,1))

#Abbildung 36: ES-Prognosen des GARCH-Modells mit der Skewed t-Verteilung
plot(dax_log_xts[1051:2050],main="ES",ylim=c(0,15))  
lines(ES995_garchsktnew_xts[1:1000],col="red")   
lines(ES99_garchsktnew_xts[1:1000],col="blue")    
lines(ES95_garchsktnew_xts[1:1000],col="green") 
legend("topleft",inset=0.005,c("ES0.005","ES0.01","ES0.05"),col=c("red","blue","green"),lty=c(1,1,1))
