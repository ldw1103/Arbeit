# Packages 
#install.packages("extRemes")
#install.packages("xts")
library(extRemes) #Extremwerttheorie
library(xts)      #eXtensible Time Series
library(tseries)
library(MASS)
library(rugarch)
library(skewt)
library(fGarch)
library(forecast)

#Daten Einlesen

daten=read.csv("DAX.csv") #Quelle: finance.yahoo.com
daten_omit=daten[!daten$Open=="null",] #Fehlende Werte wegnehmen
dim(daten_omit) #noch 6827 Beobachtungen
dax=daten_omit[,1:2]#Nur Datum und Schlusskurs bleiben
dax$Close=as.numeric(as.character(dax$Close))
dax$Date=as.Date(dax$Date)         #Datentyp von Datum transformieren
dax_xts=xts(dax$Close,dax$Date)    #xts-Format
plot(dax_xts)

n=length(dax$Close)
logreturn <- 100*log(dax$Close[-1]/dax$Close[-n])#log Return
dax_log=data.frame(dax$Date[-1],logreturn)  #2 Spalten, Datum und Log-Return
names(dax_log)[names(dax_log)=="dax.Date..1."] <- "date"  #Name verkuerzen
dax_log$date=as.Date(dax_log$date)                      #Datentyp von Datum transformieren
dax_log_xts=xts(dax_log$logreturn,dax_log$date)         #xts-Format
head(dax_log_xts)
plot(dax_log_xts)                                     #Plot von Log_Return, aber warum chinesisch?
#plot(dax_log$date,dax_log$logreturn,type = "l")
hist(dax_log_xts,breaks=100)                       #Hist von Log-Return
stats::qqnorm(dax_log$logreturn);qqline(dax_log$logreturn)  #leptokurtisch

adf.test(dax_log_xts)   #ADF-Test: p=0.01, deshalb ist es stationaer.

dax_log_xts=-dax_log_xts  # Vorzeichen umkehren, Verlust statt Rendite
plot(dax_log_xts)

#Jaehrliche Maxima als Beispiel
jahrmax_bsp=apply.yearly(dax_log_xts,max) #die jaehliche groesste Verlust 
jahrmax_bsp
plot(jahrmax_bsp)

# maximum-likelihood fitting of the GEV distribution
fit_bsp <- fevd(as.vector(jahrmax_bsp), method = "MLE", type="GEV")
# diagnostic plots
plot(fit_bsp)
fit_bsp$results$par     #Paramter. Location=3.55 (mu), Scale=1.39 (sigma), Shape=0.195 (tau) 
#xi>0, somit ist es Frechet-Verteilung

# return levels:
rl_bsp <- return.level(fit_bsp, conf = 0.05, return.period= c(2,5,10,20,50,100))
rl_bsp

#Unterschiedliche Laenge des Blocks (Woche, Monat, Quartal, Jahr)
jahrmax=apply.yearly(dax_log_xts,max) #die jaehliche groesste Verlust 
jahrmax
plot(jahrmax)
fit_jahr <- fevd(as.vector(jahrmax), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_jahr_lm <- fevd(as.vector(jahrmax), method = "Lmoments", type="GEV")  #"Lmoments" um Parameter zu schaetzen
plot(fit_jahr)
plot(fit_jahr_lm) 
return.level(fit_jahr, conf = 0.05, return.period= c(2,5,10,20,50))
return.level(fit_jahr_lm, conf = 0.05, return.period= c(2,5,10,20,50))

plot(dax_log_xts,main="Verlust und Return Level")
abline(h=return.level(fit_jahr, conf = 0.05, return.period= c(5,10,20))[1],col="green",lty=3)
abline(h=return.level(fit_jahr, conf = 0.05, return.period= c(5,10,20))[2],col="blue",lty=2)
abline(h=return.level(fit_jahr, conf = 0.05, return.period= c(5,10,20))[3],col="red",lty=1)
legend("bottomleft",inset=0.005,title="Return Level",c("20 Jahre","10 Jahre","5 Jahre"),col=c("red","blue","green"),lty=c(1,2,3))


quartalmax=apply.quarterly(dax_log_xts,max) #die quartale groesste Verlust (Empfohlen von Mcneil Calculating Quantile Risk Measures for Financial Return
                                             #Series using Extreme Value Theory)
quartalmax
plot(quartalmax)
fit_quartal <- fevd(as.vector(quartalmax), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_quartal_lm <- fevd(as.vector(quartalmax), method = "Lmoments", type="GEV")  #"Lmoments" um Parameter zu schaetzen
plot(fit_quartal)
plot(fit_quartal_lm) 
return.level(fit_quartal, conf = 0.05, return.period= c(2,5,10,20,50))
return.level(fit_quartal_lm, conf = 0.05, return.period= c(2,5,10,20,50))

monatmax=apply.monthly(dax_log_xts,max) #die monatliche groesste Verlust 
monatmax
plot(monatmax)
fit_monatmax <- fevd(as.vector(monatmax), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_monatmax_lm <- fevd(as.vector(monatmax), method = "Lmoments", type="GEV")  #"Lmoments" um Parameter zu schaetzen
plot(fit_monatmax_lm)
plot(fit_monatmax) 
return.level(fit_monatmax, conf = 0.05, return.period= c(2,5,10,20,50))
return.level(fit_monatmax_lm, conf = 0.05, return.period= c(2,5,10,20,50))


wochemax=apply.weekly(dax_log_xts,max) #die woechentliche groesste Verlust 
wochemax
plot(wochemax)
fit_woche <- fevd(as.vector(wochemax), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
plot(fit_woche) 


# Monatlich Moving Window (Groesse=1000) Zum Beispiel: das erste Fenster
ts_bm_1=dax_log_xts[1:1000]
monatmax1=apply.monthly(ts_bm_1,max)   #die groessete monatliche Verlust
monatmax1
plot(monatmax1)
fit_monat1 <- fevd(as.vector(monatmax1), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
plot(fit_monat1)         #passt
fit_monat1$results$par  #Parameter extrahieren 

# VaRs berechnen
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
monatmax=numeric(0)
fit=numeric(0)
for (i in (1:5826)){         #es gibt (6826-1000) Vorhersagen
  monatmax=apply.monthly(dax_log_xts[i:(i+999)],max)    #die groessete monatliche Verlust
  fit <- fevd(as.vector(monatmax), method = "Lmoments", type="GEV")
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 10.4833)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 20
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 5.4917)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.5588)[1]
  }
#VaR99_bmm[2040]=return.level(fevd(as.vector(-apply.monthly(dax_log_xts[2040:(2040+999)],min)), method = "Lmoment", type="GEV"), conf = 0.05, return.period= 5.255)[1]
#denn wenn i =2040, funktioniert MLE nicht, deshalb wird es dabei durch Lmoment ersetzt
VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1001:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1001:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1001:6826]) 

plot(dax_log_xts[1001:6826])  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
sum(VaR995_bmm_xts<dax_log_xts[1001:6826]) #  60 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[1001:6826]) #  109 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[1001:6826]) #  494 Ueberschreitungen

#Quartal
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
quartalmax=numeric(0)
fit=numeric(0)
for (i in (1:5826)){         #es gibt (6826-1000) Vorhersagen
  quartalmax=apply.quarterly(dax_log_xts[i:(i+999)],max)    #die groessete quartalliche Verlust
  fit <- fevd(as.vector(quartalmax), method = "Lmoments", type="GEV")
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 3.8500)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 60
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 2.2083)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.0483)[1] 
}

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1001:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1001:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1001:6826]) ##276 Problem

plot(dax_log_xts[1001:6826])  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
sum(VaR995_bmm_xts<dax_log_xts[1001:6826]) #  86 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[1001:6826]) #  164 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[1001:6826]) #  666 Ueberschreitungen

###Problem mit Daten und apply.monthly!!!!!! Die Laenge des Intervalls ist unterschiedlich!! (Das erste und das letzte!)
## Deshalb werden die Daten aequidistant unterteilt

period.max(dax_log_xts,seq(from=1,to=6826,by=20))  #Monatlich z.B.


# Moving Window (Groesse=1000) Zum Beispiel: das erste Fenster
ts_bm_1=dax_log_xts[1:1000]
monatmax1=period.max(ts_bm_1,seq(from=20,to=1000,by=20))   #die groessete monatliche Verlust
monatmax1
plot(monatmax1)
fit_monat1 <- fevd(as.vector(monatmax1), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
plot(fit_monat1)         #passt
fit_monat1$results$par  #Parameter extrahieren 

# n=20 VaRs berechnen
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
monatmax=numeric(0)
fit=numeric(0)
for (i in (1:5826)){         #es gibt (6826-1000) Vorhersagen
  monatmax=period.max(dax_log_xts[i:(i+999)],seq(from=20,to=1000,by=20))    #die groessete monatliche Verlust
  fit <- fevd(as.vector(monatmax), method = "MLE", type="GEV")
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 10.4833)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 20
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 5.4917)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.5588)[1]
}

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1001:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1001:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1001:6826]) 

plot(dax_log_xts[1001:6826])  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
sum(VaR995_bmm_xts<dax_log_xts[1001:6826]) #  63 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[1001:6826]) #  120 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[1001:6826]) #  494 Ueberschreitungen

#n=50
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
n50max=numeric(0)
fit=numeric(0)
for (i in (1:5826)){         #es gibt (6826-1000) Vorhersagen
  n50max=period.max(dax_log_xts[i:(i+999)],seq(from=50,to=1000,by=50))    #die groessete quartalliche Verlust
  fit <- fevd(as.vector(n50max), method = "MLE", type="GEV")
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 4.5109)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 50
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 2.5317)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.0834)[1] 
}

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1001:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1001:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1001:6826]) 

plot(dax_log_xts[1001:6826])  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
sum(VaR995_bmm_xts<dax_log_xts[1001:6826]) #  84 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[1001:6826]) #  160 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[1001:6826]) #  579 Ueberschreitungen


# mit Theta (extremer Index) Embrechtschap7 P.289


#GARCH
Acf(dax_log_xts)
Pacf(dax_log_xts)
adf.test(dax_log_xts)
a1=arma(dax_log_xts[1:5271],order = c(1,1))  #mit 5271 gibt es Hessian negative-semidefinite
a2=arma(dax_log_xts,order = c(1,2))
a3=arma(dax_log_xts,order = c(2,1))
a4=arma(dax_log_xts,order = c(2,2))
a5=arma(dax_log_xts,order = c(2,3))
a6=arma(dax_log_xts,order = c(3,2))

garch11 <- ugarchspec(variance.model = list(model = "sGARCH", 
                                         garchOrder = c(1, 1), 
                                         submodel = NULL, 
                                         external.regressors = NULL, 
                                         variance.targeting = FALSE), 
                   
                   mean.model     = list(armaOrder = c(1, 1), 
                                         external.regressors = NULL, 
                                         distribution.model = "norm", 
                                         start.pars = list(), 
                                         fixed.pars = list()))

garch_11 <- ugarchfit(spec = garch11, data = dax_log_xts, solver.control = list(trace=0))
garch@fit$coef
show(garch_11)
garch_11@fit$sigma
garch_11@fit$z
str(garch_11)
plot(dax_log_xts[1:250])
lines(sigma(garch_11)[1:250],col="red")




#lcation=3.55 (mu), Scale=1.39 (sigma), Shape=0.195 (xi) 
integrand1 <- function(x) {x/sigma*(1+tao*((x-mu)/sigma))^(-1/tao-1)*exp(-(1+tao*(x-mu)/sigma)^(-1/tao))}
integrand2 = function(x){mu-sigma/tao*(1-(-log((1-x)^261))^(-tao))}
sigma=1.39
mu=3.55
tao=0.195
rl_mle <- return.level(fit_mle, conf = 0.05, return.period= c(1.078,2,5,10,20,50,100))
rl_mle
integrand2(1-0.0726)
xx=integrate(integrand2,lower=0,upper=0.01)
xx$value*100
