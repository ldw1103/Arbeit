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
library(evir)
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

# Maximum-Likelihood von GEV
fit_bsp <- fevd(as.vector(jahrmax_bsp), method = "MLE", type="GEV")
# Diagnostik Plots
plot(fit_bsp)
fit_bsp$results$par     #Paramter. Location=3.55 (mu), Scale=1.39 (sigma), Shape=0.195 (tau) 
#tau >0, somit ist es Frechet-Verteilung

# return levels: (2 Jahre, 5 Jahre ...)
rl_bsp <- return.level(fit_bsp, conf = 0.05, return.period= c(2,5,10,20,50,100))
rl_bsp

#Unterschiedliche Laenge des Blocks (Jahr, Quartal, Monat, Woche)
jahrmax=apply.yearly(dax_log_xts,max) #die jaehliche groesste Verlust 
jahrmax
plot(jahrmax)
fit_jahr <- fevd(as.vector(jahrmax), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_jahr_lm <- fevd(as.vector(jahrmax), method = "Lmoments", type="GEV")  #"Lmoments" um Parameter zu schaetzen
plot(fit_jahr)
plot(fit_jahr_lm) 
return.level(fit_jahr, conf = 0.05, return.period= c(2,5,10,20,50))
return.level(fit_jahr_lm, conf = 0.05, return.period= c(2,5,10,20,50))

#Plot von Return levels anhand Jahresmaxima
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

###Problem mit Daten und apply.monthly!!!!!! Die Laenge des Intervalls ist unterschiedlich!! (Das erste und das letzte sind manchmal zu klein.)
## Deshalb werden die Daten aequidistant unterteilt wie folgt:

period.max(dax_log_xts,seq(from=1,to=6826,by=20))  #Monatlich z.B.


# Moving Window (Groesse=1200) Zum Beispiel: das erste Fenster. Laenge=1200, damit durch 20,60,30 perfekt aufteilbar
ts_bm_1=dax_log_xts[1:1200]   
monatmax1=period.max(ts_bm_1,seq(from=20,to=1200,by=20))   #die groessete monatliche Verlust
monatmax1
plot(monatmax1)
fit_monat1 <- fevd(as.vector(monatmax1), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
plot(fit_monat1)         #passt
fit_monat1$results$par  #Parameter extrahieren 

# n=20 VaRs berechnen. n=20 ist monatlich
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
monatmax=numeric(0)
fit=numeric(0)
for (i in (1:5626)){         #es gibt (6826-1200) Vorhersagen
  monatmax=period.max(dax_log_xts[i:(i+1199)],seq(from=20,to=1200,by=20))    #die groessete monatliche Verlust
  fit <- fevd(as.vector(monatmax), method = "MLE", type="GEV")
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 10.4833)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 20
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 5.4917)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.5588)[1]
}

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1201:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1201:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1201:6826]) 

plot(dax_log_xts[1201:6826])  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
sum(VaR995_bmm_xts<dax_log_xts[1201:6826]) #  67 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[1201:6826]) #  115 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[1201:6826]) #  495 Ueberschreitungen

#n=60 n =60 ist quartallich
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
n50max=numeric(0)
fit=numeric(0)
for (i in (1:5626)){         #es gibt (6826-1200) Vorhersagen
  n60max=period.max(dax_log_xts[i:(i+1199)],seq(from=60,to=1200,by=60))    
  fit <- fevd(as.vector(n60max), method = "MLE", type="GEV")
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 4.5109)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 60.
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 2.5317)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.0834)[1] 
}

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1201:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1201:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1201:6826]) 

plot(dax_log_xts[1201:6826])  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
sum(VaR995_bmm_xts<dax_log_xts[1201:6826]) #  76 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[1201:6826]) #  145 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[1201:6826]) #  495 Ueberschreitungen


# mit Theta (Extremaler Index) Embrechts 1998 Chap8 P.S19

#n=60,n=120. N_u=15,20,25,30,40,50,100,200

#n=60 
n60max=period.max(dax_log_xts,seq(from=60,to=6826,by=60))
fit_60 <- fevd(as.vector(n60max), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
#plot(fit_60)
#return.level(fit_60, conf = 0.05, return.period= c(2,5,10,20,50))
N15=quantile(dax_log_xts,(6826-15)/6826) #N_u=15
N20=quantile(dax_log_xts,(6826-20)/6826) #N_u=20
N25=quantile(dax_log_xts,(6826-25)/6826) #N_u=25
N30=quantile(dax_log_xts,(6826-30)/6826) #N_u=30
N40=quantile(dax_log_xts,(6826-40)/6826) #N_u=40
N50=quantile(dax_log_xts,(6826-50)/6826) #N_u=50
N100=quantile(dax_log_xts,(6826-100)/6826) #N_u=100
N200=quantile(dax_log_xts,(6826-200)/6826) #N_u=200
K15=sum(n60max>N15);K20=sum(n60max>N20);K25=sum(n60max>N25);K30=sum(n60max>N30);K40=sum(n60max>N40);
K50=sum(n60max>N50);K100=sum(n60max>N100);K200=sum(n60max>N200)
K15;K20;K25;K30;K40;K50;K100;K200
theta=function(n,m,Ku,Nu){(log(1-Ku/m)/log(1-Nu/6826))/n}
#Theta berechnen
theta15=theta(n=60,m=114,Ku=K15,Nu=15)
theta20=theta(n=60,m=114,Ku=K20,Nu=20)
theta25=theta(n=60,m=114,Ku=K25,Nu=25)
theta30=theta(n=60,m=114,Ku=K30,Nu=30)
theta40=theta(n=60,m=114,Ku=K40,Nu=40)
theta50=theta(n=60,m=114,Ku=K50,Nu=50)
theta100=theta(n=60,m=114,Ku=K100,Nu=100)
theta200=theta(n=60,m=114,Ku=K200,Nu=200)
theta_dach=mean(c(theta15,theta20,theta25,theta30,theta40,theta50,theta100,theta200))#Mcneil 1998 S.13,14
theta_dach #0.475  

#Theta mir "evir" berechnen. n=60,30 ...
sum(n60max>5)  #17
index60=exindex(dax_log_xts,block=60) 
index60 #0.398

#n=120 Longin 2000: Block = Semester
n120max=period.max(dax_log_xts,seq(from=120,to=6826,by=120))
sum(n120max>quantile(dax_log_xts,0.95)) #Longin: 5% als Threshold. Und Block = Semester
sum(n120max>5)
index120=exindex(dax_log_xts,block=120)
index120 #0.3190

#VaR
#n=60, Moving Window = 1200
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
ES95_bmm=numeric(0)
ES99_bmm=numeric(0)
ES995_bmm=numeric(0)
n60max=numeric(0)
fit=numeric(0)
VaR60 = function(x){mu-sigma/tau*(1-(-log((1-x)^(60*0.3190)))^(-tau))}

for (i in (1:5626)){         #es gibt (6826-1200) Vorhersagen
  n60max=period.max(dax_log_xts[i:(i+1199)],seq(from=60,to=1200,by=60))    #die groessete quartalliche Verlust
  fit <- fevd(as.vector(n60max), method = "MLE", type="GEV")
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 10.93)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^(n*theta)=(1-1/k). Hier n = 60, Theta=0.3190
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 5.7145)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.5991)[1] 
  mu=as.numeric(fit$results$par[1])
  sigma=as.numeric(fit$results$par[2])
  tau=as.numeric(fit$results$par[3])
  ES995_bmm[i]=integrate(VaR60,lower=0,upper=0.005)$value*200
  ES99_bmm[i]=integrate(VaR60,lower=0,upper=0.01)$value*100
  ES95_bmm[i]=integrate(VaR60,lower=0,upper=0.05)$value*20
}

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1201:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1201:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1201:6826]) 
ES995_bmm_xts=xts(ES995_bmm,dax_log$date[1201:6826])
ES99_bmm_xts=xts(ES99_bmm,dax_log$date[1201:6826])
ES95_bmm_xts=xts(ES95_bmm,dax_log$date[1201:6826])

plot(dax_log_xts[1201:6826],main="VaR")  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")  
legend("bottomleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))


plot(dax_log_xts[1201:6826],main="ES")  
lines(ES995_bmm_xts,col="red")   
lines(ES99_bmm_xts,col="blue")    
lines(ES95_bmm_xts,col="green")  
legend("bottomleft",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_bmm_xts<dax_log_xts[1201:6826]) #  31 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[1201:6826]) #  60 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[1201:6826]) #  248 Ueberschreitungen


#n=120
#n120max=period.max(dax_log_xts,seq(from=120,to=6826,by=120))
#fit_120 <- fevd(as.vector(n120max), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
#plot(fit_120)
#return.level(fit_120, conf = 0.05, return.period= c(2,5,10,20,50))
#N15=quantile(dax_log_xts,(6826-15)/6826) #N_u=15
#N20=quantile(dax_log_xts,(6826-20)/6826) #N_u=20
#N25=quantile(dax_log_xts,(6826-25)/6826) #N_u=25
#N30=quantile(dax_log_xts,(6826-30)/6826) #N_u=30
#N40=quantile(dax_log_xts,(6826-40)/6826) #N_u=40
#N50=quantile(dax_log_xts,(6826-50)/6826) #N_u=50
#N100=quantile(dax_log_xts,(6826-100)/6826) #N_u=100
#N200=quantile(dax_log_xts,(6826-200)/6826) #N_u=200
#K15=sum(n120max>N15);K20=sum(n120max>N20);K25=sum(n120max>N25);K30=sum(n120max>N30);K40=sum(n120max>N40);
#K50=sum(n120max>N50);K100=sum(n120max>N100);K200=sum(n120max>N200)
#K15;K20;K25;K30;K40;K50;K100;K200
#theta=function(n,m,Ku,Nu){1/n*log(1-Ku/m)/log(1-Nu/(6826))}
#Theta berechnen
#theta15=theta(n=120,m=56,Ku=K15,Nu=15)
#theta20=theta(n=120,m=56,Ku=K20,Nu=20)
#theta25=theta(n=120,m=56,Ku=K25,Nu=25)
#theta30=theta(n=120,m=56,Ku=K30,Nu=30)
#theta40=theta(n=120,m=56,Ku=K40,Nu=40)
#theta50=theta(n=120,m=56,Ku=K50,Nu=50)
#theta100=theta(n=120,m=56,Ku=K100,Nu=100)
#theta200=theta(n=120,m=56,Ku=K200,Nu=200)
#theta_dach=mean(c(theta15,theta20,theta25,theta30,theta40,theta50,theta100,theta200))
#theta_dach #0.430


###############POT################
# Mean Residual Life Plot: (Mean Excess)
mrlplot(dax_log_xts, main="Mean Residual Life Plot")    #u ist vielleicht in (0,3)
evir::meplot(dax_log_xts,xlim=c(0,6),ylim=c(1,1.5),type="l")  #u ist 3.5


# mit unterschiedlichen Grenzwerten
threshrange.plot(dax_log_xts, r = c(2, 5), nint = 16)
# ismev Implementation ist schneller:
ismev::gpd.fitrange(dax_log_xts, umin=2, umax=5, nint = 16)


# MLE
pot_mle <- fevd(as.vector(dax_log_xts), method = "MLE", type="GP", threshold=3)
# Diagnostik
plot(pot_mle)
rl_mle <- return.level(pot_mle, conf = 0.05, return.period= c(2,5,10,20,50,100))


















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
garch_11@fit$coef
show(garch_11)
garch_11@fit$sigma
garch_11@fit$z
str(garch_11)
plot(dax_log_xts[1:250])
lines(sigma(garch_11)[1:250],col="red")


