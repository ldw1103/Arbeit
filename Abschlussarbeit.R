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
plot(dax)

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

#Jaehrliche Maxima
#jahrmax=apply.yearly(dax_log_xts,min) #die jaehliche groesste Verlust 
#jahrmax=-jahrmax  #Vorzeichen umkeheren, damit sie positiv sind
#jahrmax
#plot(jahrmax)

# maximum-likelihood fitting of the GEV distribution
#fit_mle <- fevd(as.vector(jahrmax), method = "MLE", type="GEV")
# diagnostic plots
#plot(fit_mle)
#fit_mle$results$par     #Paramter. Location=3.55 (mu), Scale=1.39 (sigma), Shape=0.195 (xi) 
#xi>0, somit ist es Frechet-Verteilung

# return levels:
#rl_mle <- return.level(fit_mle, conf = 0.05, return.period= c(2,5,10,20,50,100))


#Unterschiedliche Laenge des Blocks (Woche, Monat, Quartal, Jahr)
yearmin=-apply.yearly(dax_log_xts,min) #die jaehliche groesste Verlust 
yearmin
plot(yearmin)
fit_mle <- fevd(as.vector(yearmin), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_mle <- fevd(as.vector(yearmin), method = "Lmoments", type="GEV")  #"Lmoments" um Parameter zu schaetzen
plot(fit_mle) 
plot(dax_log_xts,main="Return Level")
abline(h=return.level(fit_mle, conf = 0.05, return.period= c(5,10,20))[1],col="green",lty=3)
abline(h=return.level(fit_mle, conf = 0.05, return.period= c(5,10,20))[2],col="blue",lty=2)
abline(h=return.level(fit_mle, conf = 0.05, return.period= c(5,10,20))[3],col="red",lty=1)
legend("bottomleft",inset=0.005,title="Return Level",c("20 Jahre","10 Jahre","5 Jahre"),col=c("red","blue","green"),lty=c(1,2,3))

quartermin=-apply.quarterly(dax_log_xts,min) #die quartale groesste Verlust (Empfohlen von Mcneil Calculating Quantile Risk Measures for Financial Return
                                             #Series using Extreme Value Theory)
quartermin
plot(quartermin)
fit_mle <- fevd(as.vector(quartermin), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
plot(fit_mle) 

monthmin=-apply.monthly(dax_log_xts,min) #die monatliche groesste Verlust 
monthmin
plot(monthmin)
fit_mle <- fevd(as.vector(monthmin), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
plot(fit_mle) 

weekmin=-apply.weekly(dax_log_xts,min) #die woechentliche groesste Verlust 
weekmin
plot(weekmin)
fit_mle <- fevd(as.vector(weekmin), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
plot(fit_mle) 


# Monatlich Moving Window (Groesse=1000) Zum Beispiel: das erste Fenster
ts_bm_1=dax_log_xts[1:1000]
monthmin=-apply.monthly(ts_bm_1,min)   #die groessete monatliche Verlust
monthmin
plot(monthmin)
fit_mle <- fevd(as.vector(monthmin), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
plot(fit_mle)         #passt
fit_mle$results$par  #Parameter extrahieren 

#Alle VaRs berechnen
VaR99_bmm=numeric(0)
monthmin=numeric(0)
fit_mle=numeric(0)
for (i in (1:5826)){         #es gibt (6826-1000) Vorhersagen
  monthmin=-apply.monthly(dax_log_xts[i:(i+999)],min)    #die groessete monatliche Verlust
  fit <- fevd(as.vector(monthmin), method = "Lmoment", type="GEV")
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 5.255)[1]  #VaR0.99 Siehe Longin2000, Mcneil1998
}
#VaR99_bmm[2040]=return.level(fevd(as.vector(-apply.monthly(dax_log_xts[2040:(2040+999)],min)), method = "Lmoment", type="GEV"), conf = 0.05, return.period= 5.255)[1]
#denn wenn i =2040, funktioniert MLE nicht, deshalb wird es dabei durch Lmoment ersetzt
VaR99_bmm
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1001:6826])
VaR99_bmm_xts
plot(dax_log$logreturn)
plot(dax_log_xts[1001:6826])   
lines(VaR99_bmm_xts,col="red")    
summary(VaR99_bmm_xts>dax_log_xts[1001:6826]) #  84 Ueberschreitungen

#Quartal
VaR99_bmm=numeric(0)
quartalmin=numeric(0)
fit_mle=numeric(0)
for (i in (1:5826)){         #es gibt (6826-1000) Vorhersagen
  quartalmin=-apply.quarterly(dax_log_xts[i:(i+999)],min)    #die groessete quartalliche Verlust
  fit <- fevd(as.vector(quartalmin), method = "Lmoment", type="GEV")
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 2.085)[1]  #VaR0.99 Siehe Longin2000, Mcneil1998
}
VaR99_bmm
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1001:6826])
VaR99_bmm_xts
plot(dax_log$logreturn)
plot(dax_log_xts[1001:6826])  
lines(VaR99_bmm_xts,col="red")    
summary(VaR99_bmm_xts>dax_log_xts[1001:6826]) #  135 Ueberschreitungen




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


cc=numeric(0)
for (i in (1:5)){
cc[i]=i  
}
cc

