# Packages 
#install.packages("extRemes")
#install.packages("xts")
#install.packages("devtools")
#devtools::install_github("BayerSe/esreg")
#devtools::install_github("BayerSe/esback")
#install.packages("esback")
library(extRemes) #Extremwerttheorie
library(xts)      #eXtensible Time Series
library(tseries)  #Zeitreihe
library(MASS)
library(rugarch)  #univariate GARCH
library(skewt)    #Skewed-t-Verteilung
library(fGarch)   #GARCH
library(forecast)
library(evir)   #auch Extremwerttheorie
library(esback)  #ES backtesting

#Daten Einlesen
daten=read.csv("DAX.csv") #Quelle: finance.yahoo.com
daten_omit=daten[!daten$Open=="null",] #Fehlende Werte wegnehmen
dim(daten_omit)                     #noch 6827 Beobachtungen
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
plot(dax_log_xts)                                      #Plot von Log_Return
length(dax_log_xts)#Laenge=6826
#plot(dax_log$date,dax_log$logreturn,type = "l")
hist(dax_log_xts,breaks=100)                       #Hist von Log-Return
stats::qqnorm(dax_log$logreturn);qqline(dax_log$logreturn)  #leptokurtisch

adf.test(dax_log_xts)   #ADF-Test: p=0.01, deshalb ist es stationaer.

dax_log_xts=-dax_log_xts  # Vorzeichen umkehren, Verlust statt Rendite. --Negative Returns
plot(dax_log_xts)

mean(dax_log_xts)
median(dax_log_xts)
min(dax_log_xts)
max(dax_log_xts)
sd(dax_log_xts)
skewness(dax_log_xts)
kurtosis(dax_log_xts)

#Block-Maxima-Methode

#Jaehrliche Maxima als Beispiel
jahrmax_bsp=apply.yearly(dax_log_xts,max) #die jaehrlichen groessten Verluste 
length(jahrmax_bsp) #28 Jahre = 28 Beobachtungen
plot(jahrmax_bsp)

# Maximum-Likelihood von GEV (extRemes Package)
fit_bsp <- fevd(as.vector(jahrmax_bsp), method = "MLE", type="GEV")
# Diagnostik Plots
plot(fit_bsp)
fit_bsp$results$par     #Paramter. Location=3.55 (mu), Scale=1.39 (sigma), Shape=0.195 (xi) #xi >0, somit ist es Frechet-Verteilung
# oder mit evir Package:
fit_bsp_evir=evir::gev(dax_log_xts,block = 240)  #1 Jahre hat 240 Beobachtungen
plot(fit_bsp_evir) #gut gefittet
fit_bsp_evir$par.ests   #mu=3.44, sigma=1.44,xi=0.16


# return levels: (2 Jahre, 5 Jahre ...)
rl_bsp <- return.level(fit_bsp, conf = 0.05, return.period= c(2,5,10,20,50,100))
rl_bsp

#Unterschiedliche Laenge des Blocks (Jahr, Quartal, Monat)
jahrmax=apply.yearly(dax_log_xts,max) #die jaehliche groesste Verlust 
jahrmax
plot(jahrmax)
fit_jahr <- fevd(as.vector(jahrmax), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_jahr_lm <- fevd(as.vector(jahrmax), method = "Lmoments", type="GEV")  #"Lmoments" um Parameter zu schaetzen
plot(fit_jahr)
plot(fit_jahr_lm) 
return.level(fit_jahr, conf = 0.05, return.period= c(2,5,10,20,50))
return.level(fit_jahr_lm, conf = 0.05, return.period= c(2,5,10,20,50))

plot(dax_log_xts,main="Verlust und Return Level") #Plot von Return levels anhand Jahresmaxima
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
fit_monatmax$results$par

# Monatlich Moving Window (Groesse=1000) Zum Beispiel: das erste Fenster
#ts_bm_1=dax_log_xts[1:1000]
#monatmax1=apply.monthly(ts_bm_1,max)   #die groessete monatliche Verlust
#monatmax1
#plot(monatmax1)
#fit_monat1 <- fevd(as.vector(monatmax1), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
#plot(fit_monat1)         #passt
#fit_monat1$results$par  #Parameter extrahieren 

# VaRs berechnen
#VaR95_bmm=numeric(0)
#VaR99_bmm=numeric(0)
#VaR995_bmm=numeric(0)
#monatmax=numeric(0)
#fit=numeric(0)
#for (i in (1:5826)){         #es gibt (6826-1000) Vorhersagen
#  monatmax=apply.monthly(dax_log_xts[i:(i+999)],max)    #die groessete monatliche Verlust
#  fit <- fevd(as.vector(monatmax), method = "Lmoments", type="GEV")
#  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 10.4833)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 20
#  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 5.4917)[1]  
#  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.5588)[1]
#  }
#
#VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1001:6826])
#VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1001:6826])
#VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1001:6826]) 

#plot(dax_log_xts[1001:6826])  
#lines(VaR995_bmm_xts,col="red")   
#lines(VaR99_bmm_xts,col="blue")    
#lines(VaR95_bmm_xts,col="green")   
#sum(VaR995_bmm_xts<dax_log_xts[1001:6826]) #  60 Ueberschreitungen
#sum(VaR99_bmm_xts<dax_log_xts[1001:6826]) #  109 Ueberschreitungen
#sum(VaR95_bmm_xts<dax_log_xts[1001:6826]) #  494 Ueberschreitungen

#Quartal
#VaR95_bmm=numeric(0)
#VaR99_bmm=numeric(0)
#VaR995_bmm=numeric(0)
#quartalmax=numeric(0)
#fit=numeric(0)
#for (i in (1:5826)){         #es gibt (6826-1000) Vorhersagen
#  quartalmax=apply.quarterly(dax_log_xts[i:(i+999)],max)    #die groessete quartalliche Verlust
#  fit <- fevd(as.vector(quartalmax), method = "Lmoments", type="GEV")
#  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 3.8500)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 60
#  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 2.2083)[1]  
#  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.0483)[1] 
#}

#VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1001:6826])
#VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1001:6826])
#VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1001:6826]) ##276 Problem
#plot(dax_log_xts[1001:6826])  
#lines(VaR995_bmm_xts,col="red")   
#lines(VaR99_bmm_xts,col="blue")    
#lines(VaR95_bmm_xts,col="green")   
#sum(VaR995_bmm_xts<dax_log_xts[1001:6826]) #  86 Ueberschreitungen
#sum(VaR99_bmm_xts<dax_log_xts[1001:6826]) #  164 Ueberschreitungen
#sum(VaR95_bmm_xts<dax_log_xts[1001:6826]) #  666 Ueberschreitungen

###Problem mit z.B. Funktion "apply.monthly"!!!!!! Die Laenge des Intervalls ist unterschiedlich!! 
#(Das erste und das letzte sind manchmal zu klein!)
## Deshalb werden die Daten aequidistant unterteilt wie folgt:

period.max(dax_log_xts,seq(from=20,to=6826,by=20))  #Monatlich z.B.

###Goodness of Fit fuer GEV (Stephens 1977)
bmm_monat_st=period.max(dax_log_xts[1:6720],seq(from=20,to=6720,by=20))
bmm_quartal_st=period.max(dax_log_xts[1:6720],seq(from=60,to=6720,by=60))
bmm_semester_st=period.max(dax_log_xts[1:6720],seq(from=120,to=6720,by=120))
bmm_jahr_st=period.max(dax_log_xts[1:6720],seq(from=240,to=6720,by=240))



####Monat:
fit_monat_st <- fevd(as.vector(bmm_monat_st), method = "MLE", type="GEV") 
x_monat_st=sort(fit_monat_st$x)    ####von klein zu gross geordnet
z_monat_st=pgev(x_monat_st,xi=fit_monat_st$results$par[3],mu=fit_monat_st$results$par[1],sigma=fit_monat_st$results$par[2])

###Goodness of Fit Shermann 1957
monat_order=sort(as.numeric(bmm_monat_st))
length(monat_order) #336
pevd(monat_order[1],loc=fit_monat_st$results$par[1],scale = fit_monat_st$results$par[2],shape = fit_monat_st$results$par[3],type = "GEV")
mm=rep(0,(length(monat_order)-1))
for (i in(1:(length(monat_order)-1))){
  mm[i]=abs(pevd(monat_order[i+1],loc=fit_monat_st$results$par[1],scale = fit_monat_st$results$par[2],shape = fit_monat_st$results$par[3],type = "GEV")-
              pevd(monat_order[i],loc=fit_monat_st$results$par[1],scale = fit_monat_st$results$par[2],shape = fit_monat_st$results$par[3],type = "GEV")-
              1/(1+length(monat_order)))
}
omega=0.5*(sum(mm))+
  0.5*(abs(pevd(monat_order[1],loc=fit_monat_st$results$par[1],scale = fit_monat_st$results$par[2],shape = fit_monat_st$results$par[3],type = "GEV")-0))+
  0.5*(abs(1-pevd(monat_order[length(monat_order)],loc=fit_monat_st$results$par[1],scale = fit_monat_st$results$par[2],shape = fit_monat_st$results$par[3],type = "GEV")))
omega #statistik = 0.3676
mean_monat=(336/337)^(336+1)
var_monat=(2*exp(1)-5)/(exp(2)*336)
qnorm(0.95,mean=mean_monat,sd=sqrt(var_monat)) #Kritischer Wert= 0.3891 > 0.3676 H0 nicht abgelehnt
pnorm(omega,mean=mean_monat,sd=sqrt(var_monat)) #p= 0.51

#omega1 oder 2
omega1=1/sqrt(var_monat)*(omega-mean_monat)
omega2=sqrt(exp(2)/(2*exp(1)-1)*336)*(omega-exp(-1))
omega1 #0.020
omega2 #0.00

###Quartal
fit_quartal_st <- fevd(as.vector(bmm_quartal_st), method = "MLE", type="GEV") 
x_quartal_st=sort(fit_quartal_st$x)    ####von klein zu gross geordnet
z_quartal_st=pgev(x_quartal_st,xi=fit_quartal_st$results$par[3],mu=fit_quartal_st$results$par[1],sigma=fit_quartal_st$results$par[2])
quartal_order=sort(as.numeric(bmm_quartal_st))
length(quartal_order) #112
pevd(quartal_order[1],loc=fit_quartal_st$results$par[1],scale = fit_quartal_st$results$par[2],shape = fit_quartal_st$results$par[3],type = "GEV")
mm=rep(0,(length(quartal_order)-1))
for (i in(1:(length(quartal_order)-1))){
  mm[i]=abs(pevd(quartal_order[i+1],loc=fit_quartal_st$results$par[1],scale = fit_quartal_st$results$par[2],shape = fit_quartal_st$results$par[3],type = "GEV")-
              pevd(quartal_order[i],loc=fit_quartal_st$results$par[1],scale = fit_quartal_st$results$par[2],shape = fit_quartal_st$results$par[3],type = "GEV")-
              1/(1+length(quartal_order)))
}
omega=0.5*(sum(mm))+
  0.5*(abs(pevd(quartal_order[1],loc=fit_quartal_st$results$par[1],scale = fit_quartal_st$results$par[2],shape = fit_quartal_st$results$par[3],type = "GEV")-0))+
  0.5*(abs(1-pevd(quartal_order[length(quartal_order)],loc=fit_quartal_st$results$par[1],scale = fit_quartal_st$results$par[2],shape = fit_quartal_st$results$par[3],type = "GEV")))
omega #statistik = 0.3784
mean_quartal=(112/113)^(112+1)
var_quartal=(2*exp(1)-5)/(exp(2)*112)
qnorm(0.95,mean=mean_quartal,sd=sqrt(var_quartal)) #Kritischer Wert= 0.4040> 0.3784 H0 nicht abgelehnt
pnorm(0.3784,mean=mean_quartal,sd=sqrt(var_quartal))  #p=0.70
#omega1 oder 2
omega1=1/sqrt(var_quartal)*(omega-mean_quartal)
omega2=sqrt(exp(2)/(2*exp(1)-1)*112)*(omega-exp(-1))
omega1 #0.53
omega2 #0.14

#Semester
fit_semester_st <- fevd(as.vector(bmm_semester_st), method = "MLE", type="GEV") 
x_semester_st=sort(fit_semester_st$x)    ####von klein zu gross geordnet
z_semester_st=pgev(x_semester_st,xi=fit_semester_st$results$par[3],mu=fit_semester_st$results$par[1],sigma=fit_semester_st$results$par[2])
semester_order=sort(as.numeric(bmm_semester_st))
length(semester_order) #56
pevd(semester_order[1],loc=fit_semester_st$results$par[1],scale = fit_semester_st$results$par[2],shape = fit_semester_st$results$par[3],type = "GEV")
mm=rep(0,(length(semester_order)-1))
for (i in(1:(length(semester_order)-1))){
  mm[i]=abs(pevd(semester_order[i+1],loc=fit_semester_st$results$par[1],scale = fit_semester_st$results$par[2],shape = fit_semester_st$results$par[3],type = "GEV")-
              pevd(semester_order[i],loc=fit_semester_st$results$par[1],scale = fit_semester_st$results$par[2],shape = fit_semester_st$results$par[3],type = "GEV")-
              1/(1+length(semester_order)))
}
omega=0.5*(sum(mm))+
  0.5*(abs(pevd(semester_order[1],loc=fit_semester_st$results$par[1],scale = fit_semester_st$results$par[2],shape = fit_semester_st$results$par[3],type = "GEV")-0))+
  0.5*(abs(1-pevd(semester_order[length(semester_order)],loc=fit_semester_st$results$par[1],scale = fit_semester_st$results$par[2],shape = fit_semester_st$results$par[3],type = "GEV")))
omega #statistik = 0.3850
mean_semester=(56/57)^(56+1)
var_semester=(2*exp(1)-5)/(exp(2)*56)
qnorm(0.95,mean=mean_semester,sd=sqrt(var_semester)) #Kritischer Wert= 4180>0.3850 H0 nicht ab
pnorm(0.3850,mean=mean_semester,sd=sqrt(var_semester)) #0.7347
#omega1 oder 2
omega1=1/sqrt(var_semester)*(omega-mean_semester)
omega2=sqrt(exp(2)/(2*exp(1)-1)*56)*(omega-exp(-1))
omega1  #0.63
omega2  #0.16

##Jahr
fit_jahr_st <- fevd(as.vector(bmm_jahr_st), method = "MLE", type="GEV") 
x_jahr_st=sort(fit_jahr_st$x)    ####von klein zu gross geordnet
z_jahr_st=pgev(x_jahr_st,xi=fit_jahr_st$results$par[3],mu=fit_jahr_st$results$par[1],sigma=fit_jahr_st$results$par[2])
jahr_order=sort(as.numeric(bmm_jahr_st))
length(jahr_order) #28
pevd(jahr_order[1],loc=fit_jahr_st$results$par[1],scale = fit_jahr_st$results$par[2],shape = fit_jahr_st$results$par[3],type = "GEV")
mm=rep(0,(length(jahr_order)-1))
for (i in(1:(length(jahr_order)-1))){
  mm[i]=abs(pevd(jahr_order[i+1],loc=fit_jahr_st$results$par[1],scale = fit_jahr_st$results$par[2],shape = fit_jahr_st$results$par[3],type = "GEV")-
              pevd(jahr_order[i],loc=fit_jahr_st$results$par[1],scale = fit_jahr_st$results$par[2],shape = fit_jahr_st$results$par[3],type = "GEV")-
              1/(1+length(jahr_order)))
}
omega=0.5*(sum(mm))+
  0.5*(abs(pevd(jahr_order[1],loc=fit_jahr_st$results$par[1],scale = fit_jahr_st$results$par[2],shape = fit_jahr_st$results$par[3],type = "GEV")-0))+
  0.5*(abs(1-pevd(jahr_order[length(jahr_order)],loc=fit_jahr_st$results$par[1],scale = fit_jahr_st$results$par[2],shape = fit_jahr_st$results$par[3],type = "GEV")))
omega #statistik = 0.3244
mean_jahr=(28/29)^(28+1)
var_jahr=(2*exp(1)-5)/(exp(2)*28)
qnorm(0.95,mean=mean_jahr,sd=sqrt(var_jahr)) #Kritischer Wert= 0.4370 > 0.3244 H0 nicht ab.
pnorm(0.3244,mean=mean_jahr,sd=sqrt(var_jahr)) #0.21
#omega1 oder 2
omega1=1/sqrt(var_jahr)*(omega-mean_jahr)
omega2=sqrt(exp(2)/(2*exp(1)-1)*28)*(omega-exp(-1))
omega1 #-0.81
omega2 #-0.30

####Alle bestehen bei dem Test. Deshalb wird n = 20 gewaehlt! (mehr Beobachtungen) 



# Moving Window (Groesse=2400) Zum Beispiel: das erste Fenster. Laenge=2400, damit durch 20,60,120,240 perfekt teilbar.
ts_bm_1=dax_log_xts[1:2400]   
monatmax1=period.max(ts_bm_1,seq(from=20,to=2400,by=20))   #die groesseten monatlichen Verluste
monatmax1
plot(monatmax1)
fit_monat1 <- fevd(as.vector(monatmax1), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_monat1_evir=gev(dax_log_xts,block = 20)
plot(fit_monat1_evir)         #passt
fit_monat1$results$par  #Parameter extrahieren 

# VaRs und ESs berechnen. n=20 ist monatlich
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
ES95_bmm=numeric(0)
ES99_bmm=numeric(0)
ES995_bmm=numeric(0)
VaR20 = function(x){mu-sigma/tau*(1-(-log((1-x)^(20)))^(-tau))} #tau = xi = shape Parameter
fit=numeric(0)
#returnlevel = function(x,mu,sigma,xi){mu-sigma/xi*(1-(-log((1-x)^20))^(-xi))}
for (i in (1:4426)){         #es gibt (6826-2400) Vorhersagen
  monatmax=period.max(dax_log_xts[i:(i+2399)],seq(from=20,to=2400,by=20))    #die groesseten monatlichen Verluste
  fit <- fevd(as.vector(monatmax), method = "MLE", type="GEV")
  #fit = gev(dax_log_xts[i:(i+2399)],block=120)
  #VaR995_bmm[i]=returnlevel(0.005,mu=fit$par.ests[3],sigma=fit$par.ests[2],xi=fit$par.ests[1])#rlevel.gev(fit,2.212322)[2]
   # VaR99_bmm[i]=returnlevel(0.01,mu=fit$par.ests[3],sigma=fit$par.ests[2],xi=fit$par.ests[1])#rlevel.gev(fit,1.427308)[2]
  #  VaR95_bmm[i]=returnlevel(0.05,mu=fit$par.ests[3],sigma=fit$par.ests[2],xi=fit$par.ests[1])#rlevel.gev(fit,1.002127)[2]
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 10.483332)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 20
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 5.491697)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.558812)[1]
   mu=as.numeric(fit$results$par[1])
    sigma=as.numeric(fit$results$par[2])
    tau=as.numeric(fit$results$par[3])
    ES995_bmm[i]=integrate(VaR20,lower=0,upper=0.005)$value*200
    ES99_bmm[i]=integrate(VaR20,lower=0,upper=0.01)$value*100
    ES95_bmm[i]=integrate(VaR20,lower=0,upper=0.05)$value*20
    }

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[2401:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[2401:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[2401:6826]) 
ES995_bmm_xts=xts(ES995_bmm,dax_log$date[2401:6826])
ES99_bmm_xts=xts(ES99_bmm,dax_log$date[2401:6826])
ES95_bmm_xts=xts(ES95_bmm,dax_log$date[2401:6826])

plot(dax_log_xts[2401:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))


plot(dax_log_xts[2401:6826],main="ES",,ylim=c(0,10))  
lines(ES995_bmm_xts,col="red")   
lines(ES99_bmm_xts,col="blue")    
lines(ES95_bmm_xts,col="green")  
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_bmm_xts<dax_log_xts[2401:6826]) # 50 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[2401:6826]) #  83 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[2401:6826]) #  387 Ueberschreitungen

#U.C Test (Code aus dem Buch von Danielsson)
uc_test = function(p,v){
  return(2*log(((1-sum(v)/length(v))/(1-p))^(length(v)-sum(v))*((sum(v)/length(v))/p)^(sum(v))))
}

V995=(VaR995_bmm_xts<dax_log_xts[2401:6826])
sum(V995)
V99=(VaR99_bmm_xts<dax_log_xts[2401:6826])
sum(V99)
V95=(VaR95_bmm_xts<dax_log_xts[2401:6826])
sum(V95)

uc_test(p=0.005,v=V995) #25.94
uc_test(p=0.01,v=V99)  #27.24
uc_test(p=0.05,v=V95)  #107.8 #2*log(((1-387/4426)/0.95)^(4426-387)*((387/4426)/0.05)^(387))=107.8
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)
1-pchisq(uc_test(p=0.01,v=V99),1)
1-pchisq(uc_test(p=0.05,v=V95),1)

##Ind.Test  (Code aus dem Buch von Danielsson)
ind_test = function(V){
  J = matrix(ncol = 4,nrow = length(V))
  for (i in 2:length(V)){
    J[i,1] = V[i - 1] == 0 & V[i] == 0
    J[i,2] = V[i - 1] ==0 & V[i] == 1
    J[i,3] = V[i - 1] == 1 & V[i] == 0
    J[i,4] = V[i - 1] == 1 & V[i] == 1
  }
  V_00 = sum(J[,1],na.rm = TRUE)
  V_01 = sum(J[,2],na.rm = TRUE)
  V_10 = sum(J[,3],na.rm = TRUE)
  V_11 = sum(J[,4],na.rm = TRUE)
  p_00 = V_00/(V_00 + V_01)
  p_01 = V_01/(V_00 + V_01)
  p_10 = V_10/(V_10 + V_11)
  p_11 = V_11/(V_10 + V_11)
  hat_p = (V_01 + V_11)/(V_00 + V_01 + V_10 + V_11)
  #a = (1 - hat_p)^(V_00 + V_10)*(hat_p)^(V_01 + V_11)
  #b = (p_00)^(V_00)*(p_01)^(V_01)*(p_10)^(V_10)* p_11^(V_11)
  return(-2 * log(((1 - hat_p)/p_00)^(V_00)*((1 - hat_p)/p_10)^(V_10)*(hat_p/p_01)^(V_01)*(hat_p/p_11)^(V_11)))
}

ind_test(as.vector(V995))#5.39
ind_test(as.vector(V99)) #14.36
ind_test(as.vector(V95))  #41.34
1-pchisq(ind_test(as.vector(V995)),1)#0.02
1-pchisq(ind_test(as.vector(V99)),1)#0.00
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.00
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.00
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00

##ES Test (Mcneil 2000)
#ESTest(0.005,-dax_log_xts[2401:6826],ES=-ES995_bmm_xts,VaR=-VaR995_bmm_xts)
#ESTest(0.01,-dax_log_xts[2401:6826],ES=-ES99_bmm_xts,VaR=-VaR99_bmm_xts)
#ESTest(0.05,actual=-dax_log_xts[2401:6826],ES=-ES95_bmm_xts,VaR=-VaR95_bmm_xts)

####Acerbi Test 2
Z2=function(p,ES,L,v){
  s = matrix(ncol = 1,nrow = length(ES))
  for (i in 1:length(ES)){
  s[i]=L[i]*v[i]/(p*length(ES)*ES[i])
  }
  return(sum(s)-1)
}
Z2(p=0.005,ES=ES995_bmm_xts,L=dax_log_xts[2401:6826],v=V995)#1.04
Z2(p=0.01,ES=ES99_bmm_xts,L=dax_log_xts[2401:6826],v=V99)#0.84
Z2(p=0.05,ES=ES95_bmm_xts,L=dax_log_xts[2401:6826],v=V95)#0.77

#ESBACK
#esback::esr_backtest(-dax_log_xts[2401:6826],e=-ES995_bmm_xts,alpha=0.005)
#esback::esr_backtest(-dax_log_xts[2401:6826],e=-ES99_bmm_xts,alpha=0.01)
#esback::esr_backtest(-dax_log_xts[2401:6826],e=-ES95_bmm_xts,alpha=0.05)

esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES995_bmm_xts,alpha=0.005)#0.12
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES99_bmm_xts,alpha=0.01)#0.00
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES95_bmm_xts,alpha=0.05)#0.00

#esback::cc_backtest(-dax_log_xts[2401:6826],q=-ES995_bmm_xts,e=-VaR995_bmm_xts,alpha=0.005)
#esback::cc_backtest(-dax_log_xts[2401:6826],q=-ES99_bmm_xts,e=-VaR99_bmm_xts,alpha=0.01)
#esback::cc_backtest(-dax_log_xts[2401:6826],q=-ES99_bmm_xts,e=-VaR95_bmm_xts,alpha=0.05)

#esback::er_backtest(-dax_log_xts[2401:6826],q=-ES995_bmm_xts,e=-VaR995_bmm_xts)
#esback::er_backtest(-dax_log_xts[2401:6826],q=-ES99_bmm_xts,e=-VaR99_bmm_xts)
#esback::er_backtest(-dax_log_xts[2401:6826],q=-ES95_bmm_xts,e=-VaR95_bmm_xts)



# Mit Theta (Extremaler Index) Embrechts 1998 Chap8 P.S19

#n=60,n=120. N_u=15,20,25,30,40,50,100,200

#n=60  #Mcneil 1998 S.13,14  (m, n muessen gross genug sein)
n60max=period.max(dax_log_xts[1:6720],seq(from=60,to=6720,by=60))
fit_60 <- fevd(as.vector(n60max), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
#plot(fit_60)
#return.level(fit_60, conf = 0.05, return.period= c(2,5,10,20,50))
N15=quantile(dax_log_xts,(6720-15)/6720) #N_u=15
N20=quantile(dax_log_xts,(6720-20)/6720) #N_u=20
N25=quantile(dax_log_xts,(6720-25)/6720) #N_u=25
N30=quantile(dax_log_xts,(6720-30)/6720) #N_u=30
N40=quantile(dax_log_xts,(6720-40)/6720) #N_u=40
N50=quantile(dax_log_xts,(6720-50)/6720) #N_u=50
N100=quantile(dax_log_xts,(6720-100)/6720) #N_u=100
N200=quantile(dax_log_xts,(6720-200)/6720) #N_u=200
K15=sum(n60max>N15);K20=sum(n60max>N20);K25=sum(n60max>N25);K30=sum(n60max>N30);K40=sum(n60max>N40);
K50=sum(n60max>N50);K100=sum(n60max>N100);K200=sum(n60max>N200)
K15;K20;K25;K30;K40;K50;K100;K200
theta=function(n,m,Ku,Nu){(log(1-Ku/m)/log(1-Nu/6720))/n}
#Theta berechnen
theta15=theta(n=60,m=112,Ku=K15,Nu=15)
theta20=theta(n=60,m=112,Ku=K20,Nu=20)
theta25=theta(n=60,m=112,Ku=K25,Nu=25)
theta30=theta(n=60,m=112,Ku=K30,Nu=30)
theta40=theta(n=60,m=112,Ku=K40,Nu=40)
theta50=theta(n=60,m=112,Ku=K50,Nu=50)
theta100=theta(n=60,m=112,Ku=K100,Nu=100)
theta200=theta(n=60,m=112,Ku=K200,Nu=200)
theta_dach=mean(c(theta15,theta20,theta25,theta30,theta40,theta50,theta100,theta200))#Mcneil 1998 S.13,14
theta_dach #0.50  

#n=120
n120max=period.max(dax_log_xts[1:6720],seq(from=120,to=6720,by=120))
fit_120 <- fevd(as.vector(n120max), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
#plot(fit_60)
#return.level(fit_60, conf = 0.05, return.period= c(2,5,10,20,50))
N15=quantile(dax_log_xts,(6720-15)/6720) #N_u=15
N20=quantile(dax_log_xts,(6720-20)/6720) #N_u=20
N25=quantile(dax_log_xts,(6720-25)/6720) #N_u=25
N30=quantile(dax_log_xts,(6720-30)/6720) #N_u=30
N40=quantile(dax_log_xts,(6720-40)/6720) #N_u=40
N50=quantile(dax_log_xts,(6720-50)/6720) #N_u=50
N100=quantile(dax_log_xts,(6720-100)/6720) #N_u=100
N200=quantile(dax_log_xts,(6720-200)/6720) #N_u=200
K15=sum(n120max>N15);K20=sum(n120max>N20);K25=sum(n120max>N25);K30=sum(n120max>N30);K40=sum(n120max>N40);
K50=sum(n120max>N50);K100=sum(n120max>N100);K200=sum(n120max>N200)
K15;K20;K25;K30;K40;K50;K100;K200
theta=function(n,m,Ku,Nu){(log(1-Ku/m)/log(1-Nu/6720))/n}
#Theta berechnen
theta15=theta(n=120,m=56,Ku=K15,Nu=15)
theta20=theta(n=120,m=56,Ku=K20,Nu=20)
theta25=theta(n=120,m=56,Ku=K25,Nu=25)
theta30=theta(n=120,m=56,Ku=K30,Nu=30)
theta40=theta(n=120,m=56,Ku=K40,Nu=40)
theta50=theta(n=120,m=56,Ku=K50,Nu=50)
theta100=theta(n=120,m=56,Ku=K100,Nu=100)
theta200=theta(n=120,m=56,Ku=K200,Nu=200)
theta_dach=mean(c(theta15,theta20,theta25,theta30,theta40,theta50,theta100,theta200))#Mcneil 1998 S.13,14
theta_dach #0.44  



#Theta mir "evir" berechnen. n=60,120 ...

#n=120 Longin 2000: Block = Semester
#n120max=period.max(dax_log_xts[1:6720],seq(from=120,to=6720,by=120))
#sum(n120max>quantile(dax_log_xts,0.95))
#sum(n120max>5) #Ueberschreitungen des Blocks = 13. Longin 2000: 5% als Threshold. Und Block = Semester
#index120=exindex(dax_log_xts,block=120)
#index120 #0.3190

#sum(n60max>5)  #17  #Longin 2000: 5% als Threshold. 
#index60=exindex(dax_log_xts,block=60) 
#index60 #0.398




#VaR
#n=120, Moving Window = 2400. Theta= 0.50. (S.13 Mcneil1998)
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
ES95_bmm=numeric(0)
ES99_bmm=numeric(0)
ES995_bmm=numeric(0)
n20max=numeric(0)
fit=numeric(0)
VaR20 = function(x){mu-sigma/tau*(1-(-log((1-x)^(20*0.50)))^(-tau))} #tau = xi = shape Parameter

for (i in (1:4426)){         #es gibt (6826-2400) Vorhersagen
  n20max=period.max(dax_log_xts[i:(i+2399)],seq(from=20,to=2400,by=20))    #die groessete quartalliche Verlust
  fit <- fevd(as.vector(n20max), method = "MLE", type="GEV")
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 20.454135)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^(n*theta)=(1-1/k). Hier n = 20, Theta=0.45
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 10.458290)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 2.492131)[1] 
  mu=as.numeric(fit$results$par[1])
  sigma=as.numeric(fit$results$par[2])
  tau=as.numeric(fit$results$par[3])
  ES995_bmm[i]=integrate(VaR20,lower=0,upper=0.005)$value*200
  ES99_bmm[i]=integrate(VaR20,lower=0,upper=0.01)$value*100
  ES95_bmm[i]=integrate(VaR20,lower=0,upper=0.05)$value*20
}

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[2401:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[2401:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[2401:6826]) 
ES995_bmm_xts=xts(ES995_bmm,dax_log$date[2401:6826])
ES99_bmm_xts=xts(ES99_bmm,dax_log$date[2401:6826])
ES95_bmm_xts=xts(ES95_bmm,dax_log$date[2401:6826])

plot(dax_log_xts[2401:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")  
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))


plot(dax_log_xts[2401:6826],main="ES",ylim=c(0,10))  
lines(ES995_bmm_xts,col="red")   
lines(ES99_bmm_xts,col="blue")    
lines(ES95_bmm_xts,col="green")  
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_bmm_xts<dax_log_xts[2401:6826]) #  20 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[2401:6826]) #  50 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[2401:6826]) #  204 Ueberschreitungen

#U.C Test
V995=(VaR995_bmm_xts<dax_log_xts[2401:6826])
sum(V995)
V99=(VaR99_bmm_xts<dax_log_xts[2401:6826])
sum(V99)
V95=(VaR95_bmm_xts<dax_log_xts[2401:6826])
sum(V95)

uc_test(p=0.005,v=V995) #0.213
uc_test(p=0.01,v=V99)  #0.722
uc_test(p=0.05,v=V95)  #1.46
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)
1-pchisq(uc_test(p=0.01,v=V99),1)
1-pchisq(uc_test(p=0.05,v=V95),1)

##Ind.Test 
ind_test(as.vector(V995))#3.07
ind_test(as.vector(V99)) #5.39
ind_test(as.vector(V95))  #35.86
1-pchisq(ind_test(as.vector(V995)),1)#0.08
1-pchisq(ind_test(as.vector(V99)),1)#0.02
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.19
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.05
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00

##ES Test (Mcneil 2000)
#ESTest(0.005,-dax_log_xts[2401:6826],ES=-ES995_bmm_xts,VaR=-VaR995_bmm_xts)
#ESTest(0.01,-dax_log_xts[2401:6826],ES=-ES99_bmm_xts,VaR=-VaR99_bmm_xts)
#ESTest(0.05,actual=-dax_log_xts[2401:6826],ES=-ES95_bmm_xts,VaR=-VaR95_bmm_xts)

####Acerbi Test 2
Z2=function(p,ES,L,v){
  s = matrix(ncol = 1,nrow = length(ES))
  for (i in 1:length(ES)){
    s[i]=L[i]*v[i]/(p*length(ES)*ES[i])
  }
  return(sum(s)-1)
}
Z2(p=0.005,ES=ES995_bmm_xts,L=dax_log_xts[2401:6826],v=V995)#-0.17
Z2(p=0.01,ES=ES99_bmm_xts,L=dax_log_xts[2401:6826],v=V99)#0.02
Z2(p=0.05,ES=ES95_bmm_xts,L=dax_log_xts[2401:6826],v=V95)#-0.08

#ESBACK
#esback::esr_backtest(-dax_log_xts[2401:6826],e=-ES995_bmm_xts,alpha=0.005)
#esback::esr_backtest(-dax_log_xts[2401:6826],e=-ES99_bmm_xts,alpha=0.01)
#esback::esr_backtest(-dax_log_xts[2401:6826],e=-ES95_bmm_xts,alpha=0.05)

esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES995_bmm_xts,alpha=0.005)#0.99
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES99_bmm_xts,alpha=0.01)#0.96
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES95_bmm_xts,alpha=0.05)#0.83

#esback::cc_backtest(-dax_log_xts[2401:6826],q=-ES995_bmm_xts,e=-VaR995_bmm_xts,alpha=0.005)
#esback::cc_backtest(-dax_log_xts[2401:6826],q=-ES99_bmm_xts,e=-VaR99_bmm_xts,alpha=0.01)
#esback::cc_backtest(-dax_log_xts[2401:6826],q=-ES99_bmm_xts,e=-VaR95_bmm_xts,alpha=0.05)

#esback::er_backtest(-dax_log_xts[2401:6826],q=-ES995_bmm_xts,e=-VaR995_bmm_xts)
#esback::er_backtest(-dax_log_xts[2401:6826],q=-ES99_bmm_xts,e=-VaR99_bmm_xts)
#esback::er_backtest(-dax_log_xts[2401:6826],q=-ES95_bmm_xts,e=-VaR95_bmm_xts)





###############POT################
# Mean Residual Life Plot: (Mean Excess)
mrlplot(dax_log_xts, main="Mean Residual Life Plot",)    #u ist vielleicht in (0,4)
meplot(dax_log_xts,xlim=c(0,5),ylim=c(1,1.5),type="l")  #u ist vielleicht 3.5. Nach 3.5 ist linear. Ist (6826-100)/6826=0.9854Quantil


#Hill-Schaetzer (Hill, 1975) (Mecneil 2000) tau^dach=1/k*Sigma^k_(j=1)(log(zj)-log(z(k+1))). 

#####Aber bei Hill-Schaetzer: shape-Parameter muss >0! (Mcneil 2000 Seite.288) Kann noch als Instrument zur Threshold-Wahl?
n=6826
hill(dax_log_xts)
evir::hill(dax_log_xts,xlim=c(60,340))  #Hill Plot.Ab 100 ist es linear  k=ungefaehr 100, y-Achse = ungefaehr 3.2.
quantile(dax_log_xts,(6826-100)/6826)  # Threshold wird als 3.50 gewaehlt.


#dax_log_order=sort(-dax_log$logreturn,decreasing = TRUE) 
#dax_log_order[100]#u = 3.50

taudach=numeric(0)
for (i in (60:340)){
  taudach[i]=1/i*sum(log(dax_log_order[1:i])-log(dax_log_order[i]))
}
plot(taudach^(-1),type="l")    #identisch zur evir::hill. Ab ungefaehr 100 ise es linear

# mit unterschiedlichen Grenzwerten
threshrange.plot(dax_log_xts, r = c(0, 5), nint = 16,type="GP")
# ismev Implementation ist schneller:
ismev::gpd.fitrange(dax_log_xts, umin=0, umax=5, nint = 50) 


# MLE mit extRemes
pot_mle <- fevd(as.vector(dax_log_xts), method = "MLE", type="GP", threshold=3.50)
# Diagnostik
pot_mle$results$par  #geschaetzte Parameter
plot(pot_mle)

#Mle mit evir
pot_mle_evir=gpd(dax_log_xts,threshold=3.50,method = "ml")
pot_mle_evir$par.ests
plot(pot_mle_evir) #diagnostik. gut gepasst

#Unbedingte VaRs-Schaetzung
r=riskmeasures(pot_mle_evir,c(0.95,0.99,0.995))
r

#Moving Windows mit Laenge 2400
VaR95_pot=numeric(0)
VaR99_pot=numeric(0)
VaR995_pot=numeric(0)
ES95_pot=numeric(0)
ES99_pot=numeric(0)
ES995_pot=numeric(0)
for (i in (1:4426)){         #es gibt (6826-2400) Vorhersagen  
  gpdpot=gpd(dax_log_xts[i:(2399+i)],threshold=quantile(dax_log_xts[i:(2399+i)],0.9),method = "ml")
  #Mcneil,2000 S.288 the choice of k in Moving Window (10%)
  risk=riskmeasures(gpdpot,c(0.95,0.99,0.995))
  VaR995_pot[i]=risk[6]
  VaR99_pot[i]=risk[5]
  VaR95_pot[i]=risk[4]
  ES995_pot[i]=risk[9]
  ES99_pot[i]=risk[8]
  ES95_pot[i]=risk[7]
}  

VaR995_pot_xts=xts(VaR995_pot,dax_log$date[2401:6826])
VaR99_pot_xts=xts(VaR99_pot,dax_log$date[2401:6826])
VaR95_pot_xts=xts(VaR95_pot,dax_log$date[2401:6826]) 

ES995_pot_xts=xts(ES995_pot,dax_log$date[2401:6826])
ES99_pot_xts=xts(ES99_pot,dax_log$date[2401:6826])
ES95_pot_xts=xts(ES95_pot,dax_log$date[2401:6826])


plot(dax_log_xts[2401:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_pot_xts,col="red")   
lines(VaR99_pot_xts,col="blue")    
lines(VaR95_pot_xts,col="green") 
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

plot(dax_log_xts[2401:6826],main="ES",ylim=c(0,10))  
lines(ES995_pot_xts,col="red")   
lines(ES99_pot_xts,col="blue")    
lines(ES95_pot_xts,col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_pot_xts<dax_log_xts[2401:6826]) #  26 Ueberschreitungen
sum(VaR99_pot_xts<dax_log_xts[2401:6826]) #  61 Ueberschreitungen
sum(VaR95_pot_xts<dax_log_xts[2401:6826]) #  222 Ueberschreitungen

#U.C Test
V995=(VaR995_pot_xts<dax_log_xts[2401:6826])
sum(V995)
V99=(VaR99_pot_xts<dax_log_xts[2401:6826])
sum(V99)
V95=(VaR95_pot_xts<dax_log_xts[2401:6826])
sum(V95)

uc_test(p=0.005,v=V995) #0.6438
uc_test(p=0.01,v=V99)  #5.7201
uc_test(p=0.05,v=V95)  #0.0023
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.4223216
1-pchisq(uc_test(p=0.01,v=V99),1)#0.01676582
1-pchisq(uc_test(p=0.05,v=V95),1)#0.9615142

##Ind.Test 
ind_test(as.vector(V995))#6.864794
ind_test(as.vector(V99)) #6.498452
ind_test(as.vector(V95))  #27.70598
1-pchisq(ind_test(as.vector(V995)),1)#0.01
1-pchisq(ind_test(as.vector(V99)),1)#0.01
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.02
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.00
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00


####Acerbi Test 2
Z2=function(p,ES,L,v){
  s = matrix(ncol = 1,nrow = length(ES))
  for (i in 1:length(ES)){
    s[i]=L[i]*v[i]/(p*length(ES)*ES[i])
  }
  return(sum(s)-1)
}
Z2(p=0.005,ES=ES995_pot_xts,L=dax_log_xts[2401:6826],v=V995)#0.17
Z2(p=0.01,ES=ES99_pot_xts,L=dax_log_xts[2401:6826],v=V99)#0.34
Z2(p=0.05,ES=ES95_pot_xts,L=dax_log_xts[2401:6826],v=V95)#0.02

esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES995_pot_xts,alpha=0.005)#0.33
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES99_pot_xts,alpha=0.01)#0.24
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES95_pot_xts,alpha=0.05)#0.24


######Danielsson2001 Threshold mit Hilfe von Subsample-Bootstrap
#install.packages("tea")
#library(tea) # Package zum Danielssons Bootstrap
#danielsson(dax_log$logreturn,B=100) #Threshold= 4.20. aber das einzelne Verfahren kostet mehr als 20 Minuten

#eye(dax_log$logreturn)  #error

#GH(dax_log$logreturn)  #funktioniert. Guillou,A.andHall,P.(2001)ADiagnosticforSelectingtheThresholdinExtremeValueAnalysis

#gomes(dax_log$logreturn,B=10,epsilon = 0.995)  #error

#hall(dax_log$logreturn)#Threshold = 2.86

#Himp(dax_log$logreturn)

#Moving Windows mit Laenge 2400
#VaR95_pot=numeric(0)
#VaR99_pot=numeric(0)
#VaR995_pot=numeric(0)
#ES95_pot=numeric(0)
#ES99_pot=numeric(0)
#ES995_pot=numeric(0)
#for (i in (1:4426)){         #es gibt (6826-2400) Vorhersagen  
#  gg=GH(dax_log$logreturn[i:(2399+i)])
#  ts=gg$threshold
#  gpdpot=gpd(dax_log$logreturn[i:(2399+i)],threshold=ts,method = "ml")
#  risk=riskmeasures(gpdpot,c(0.95,0.99,0.995))
#  VaR995_pot[i]=risk[6]
#  VaR99_pot[i]=risk[5]
#  VaR95_pot[i]=risk[4]
#  ES995_pot[i]=risk[9]
#  ES99_pot[i]=risk[8]
#  ES95_pot[i]=risk[7]
#}  ###non-finite finite-difference value [1]

#VaR995_pot_xts=xts(VaR995_pot,dax_log$date[2401:6826])
#VaR99_pot_xts=xts(VaR99_pot,dax_log$date[2401:6826])
#VaR95_pot_xts=xts(VaR95_pot,dax_log$date[2401:6826]) 

#plot(dax_log_xts[2401:6826],main="VaR")  
#lines(VaR995_pot_xts,col="red")   
#lines(VaR99_pot_xts,col="blue")    
#lines(VaR95_pot_xts,col="green")   
#legend("bottomleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#sum(VaR995_pot_xts<dax_log_xts[1201:6826]) #  41 Ueberschreitungen
#sum(VaR99_pot_xts<dax_log_xts[1201:6826]) #  76 Ueberschreitungen
#sum(VaR95_pot_xts<dax_log_xts[1201:6826]) #  316 Ueberschreitungen



##POT mit GARCH-Filter

library(fGarch)
garchfit1=garchFit(formula=~garch(1,1),data=dax_log_xts)
garchfit1

garchfit2=garchFit(formula=~garch(1,1),data=dax_log_xts,cond.dist ="QMLE") #Mcneil2000 Pseudo-MLE
garchfit2

#####AR(1)-GARCH(1,1) Mcneil 2000
garchfit3=garchFit(formula=~arma(1,0)+garch(1,1),data=dax_log_xts,cond.dist ="QMLE")
garchfit3 #ar1 nicht sig. deshalb arma(0,0)
garchfit3=garchFit(formula=~garch(1,1),data=dax_log_xts,cond.dist ="QMLE")
garchfit3
garchfit3@fit$par
garchfit3@residuals
garchfit3@sigma.t    #sd
garchfit3@h.t       #Var  
garchfit3@fitted    

zt=(dax_log_xts-garchfit3@fitted)/garchfit3@sigma.t   #standardisierte Residuen
plot(zt,main="Standardisierte Residuen")
plot(abs(dax_log_xts),main="Betrag der Log-Verluste")
plot(garchfit3@sigma.t,type="l",main="Standardabweichung der Reihe")

par(mfrow=c(2,2))
acf(dax_log_xts,main="ACF der Log-Verluste");acf(abs(dax_log_xts),main="ACF der Abs(Log-Verluste)")  #nicht i.i.d
acf(zt,main="ACF der Standardisierten Residuen");acf(abs(zt),main="ACF der Abs(Standardisierten Residuen)")     #keinen ARCH-Effekt

Box.test(dax_log_xts,lag=10,type="Ljung-Box") #p=0.00 H0: unabhaengig 
Box.test(zt,lag=10,type="Ljung-Box")  #p=0.53

stats::qqnorm(zt);qqline(zt)

# Mean Residual Life Plot: (Mean Excess)
mrlplot(zt, main="Mean Residual Life Plot")    
meplot(zt,type="l")  
meplot(zt,type="l")  #Threshold = ungefaehr 2.5. Ab 2.5 ist es linear


#Hill-Plot
evir::hill(zt,xlim=c(15,340))  #Hill Plot  k=ungefaehr 70, y = ungefaehr 5.
quantile(zt,(6826-70)/6826)  # Threshold wird als 2.54 gewaehlt.


zt_order=sort((-dax_log$logreturn-garchfit3@fitted)/garchfit3@sigma.t,decreasing = TRUE) 
zt_order[70]#u = 2.54

taudach1=numeric(0)
for (i in (15:500)){
  taudach1[i]=1/i*sum(log(zt_order[1:i])-log(zt_order[i]))
}
plot(taudach1,type="l")    #identisch zur evir::hill. 

# mit unterschiedlichen Grenzwerten
threshrange.plot(zt, r = c(2, 4), nint = 16)
# ismev Implementation ist schneller:
ismev::gpd.fitrange(zt, umin=2, umax=4, nint = 16) # nicht informativ


# MLE mit extRemes
pot_mle_garch <- fevd(as.vector(zt), method = "MLE", type="GP", threshold=2.54)
# Diagnostik
pot_mle_garch$results$par  #geschaetzte Parameter
plot(pot_mle_garch)

#Mle mit evir
pot_mle_evir_garch=gpd(zt,threshold=2.54,method = "ml")
pot_mle_evir_garch$par.ests
#par(mfrow=c(2,2))
plot(pot_mle_evir_garch) #diagnostik. gut gepasst

#Unbedingte VaR-Schaetzung
r=riskmeasures(pot_mle_evir_garch,c(0.95,0.99,0.995))
r

#Moving Windows mit Laenge 2400 (Fuer jedes MovingWindow wird GARCH erneut geschaetzt)
VaR95_pot_z=numeric(0)
VaR99_pot_z=numeric(0)
VaR995_pot_z=numeric(0)
ES95_pot_garch=numeric(0)
ES99_pot_garch=numeric(0)
ES995_pot_garch=numeric(0)

VaR95_pot_garch=numeric(0)
VaR99_pot_garch=numeric(0)
VaR995_pot_garch=numeric(0)
ES95_pot_garch=numeric(0)
ES99_pot_garch=numeric(0)
ES995_pot_garch=numeric(0)

for (i in (1:4426)){         #es gibt (6826-2400) Vorhersagen. (laeuft 45 Minuten..)
  garchfitm=garchFit(formula=~arma(1,0)+garch(1,1),data=dax_log_xts[i:(2399+i)],cond.dist ="QMLE")
  ztm=(dax_log_xts[i:(2399+i)]-garchfitm@fitted)/garchfitm@sigma.t
  gpdpotgarch=fevd(as.vector(ztm), method = "MLE", type="GP", threshold=quantile(ztm,(2400-240)/2400))
  VaR995_pot_z=quantile(ztm,0.9)+gpdpotgarch$results$par[1]/gpdpotgarch$results$par[2]*((2400*0.005/240)^(-gpdpotgarch$results$par[2])-1)
  VaR99_pot_z=quantile(ztm,0.9)+gpdpotgarch$results$par[1]/gpdpotgarch$results$par[2]*((2400*0.01/240)^(-gpdpotgarch$results$par[2])-1)
  VaR95_pot_z=quantile(ztm,0.9)+gpdpotgarch$results$par[1]/gpdpotgarch$results$par[2]*((2400*0.05/240)^(-gpdpotgarch$results$par[2])-1)
  
  VaR995_pot_garch[i]=VaR995_pot_z*predict(garchfitm)[1,3]+predict(garchfitm)[1,1]#Mcneil S.6
  VaR99_pot_garch[i]=VaR99_pot_z*predict(garchfitm)[1,3]+predict(garchfitm)[1,1]
  VaR95_pot_garch[i]=VaR95_pot_z*predict(garchfitm)[1,3]+predict(garchfitm)[1,1]
  #ES: Mcneil S.293
  ES995_pot_garch[i]=predict(garchfitm)[1,1]+predict(garchfitm)[1,3]*VaR995_pot_z*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(ztm,0.9))/(VaR995_pot_z-VaR995_pot_z*gpdpotgarch$results$par[2]))
  ES99_pot_garch[i]=predict(garchfitm)[1,1]+predict(garchfitm)[1,3]*VaR99_pot_z*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(ztm,0.9))/(VaR99_pot_z-VaR995_pot_z*gpdpotgarch$results$par[2]))
  ES95_pot_garch[i]=predict(garchfitm)[1,1]+predict(garchfitm)[1,3]*VaR95_pot_z*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(ztm,0.9))/(VaR95_pot_z-VaR995_pot_z*gpdpotgarch$results$par[2]))
}  
VaR995_pot_xts_garch=xts(VaR995_pot_garch,dax_log$date[2401:6826])
VaR99_pot_xts_garch=xts(VaR99_pot_garch,dax_log$date[2401:6826])
VaR95_pot_xts_garch=xts(VaR95_pot_garch,dax_log$date[2401:6826])

ES995_pot_xts_garch=xts(ES995_pot_garch,dax_log$date[2401:6826])
ES99_pot_xts_garch=xts(ES99_pot_garch,dax_log$date[2401:6826])
ES95_pot_xts_garch=xts(ES95_pot_garch,dax_log$date[2401:6826])



plot(dax_log_xts[2401:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_pot_xts_garch,col="red")   
lines(VaR99_pot_xts_garch,col="blue")    
lines(VaR95_pot_xts_garch,col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_pot_xts_garch<dax_log_xts[2401:6826]) #  16 Ueberschreitungen
sum(VaR99_pot_xts_garch<dax_log_xts[2401:6826]) #  36 Ueberschreitungen
sum(VaR95_pot_xts_garch<dax_log_xts[2401:6826]) #  232 Ueberschreitungen

plot(dax_log_xts[2401:6826],main="ES",ylim=c(0,10))  
lines(ES995_pot_xts_garch,col="red")   
lines(ES99_pot_xts_garch,col="blue")    
lines(ES95_pot_xts_garch,col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))


###ARMA(0,0) GARCH(1,1)
VaR95_pot_z0=numeric(0)
VaR99_pot_z0=numeric(0)
VaR995_pot_z0=numeric(0)


VaR95_pot_garch0=numeric(0)
VaR99_pot_garch0=numeric(0)
VaR995_pot_garch0=numeric(0)
ES95_pot_garch0=numeric(0)
ES99_pot_garch0=numeric(0)
ES995_pot_garch0=numeric(0)


for (i in (1:4426)){         #es gibt (6826-2400) Vorhersagen. (laeuft 45 Minuten..)
  garchfitm0=garchFit(formula=~arma(0,0)+garch(1,1),data=dax_log_xts[i:(2399+i)],cond.dist ="QMLE")
  ztm0=(dax_log_xts[i:(2399+i)]-garchfitm0@fitted)/garchfitm0@sigma.t
  gpdpotgarch0=fevd(as.vector(ztm0), method = "MLE", type="GP", threshold=quantile(ztm0,(2400-240)/2400))
  VaR995_pot_z0=quantile(ztm0,0.9)+gpdpotgarch0$results$par[1]/gpdpotgarch0$results$par[2]*((2400*0.005/240)^(-gpdpotgarch0$results$par[2])-1)
  VaR99_pot_z0=quantile(ztm0,0.9)+gpdpotgarch0$results$par[1]/gpdpotgarch0$results$par[2]*((2400*0.01/240)^(-gpdpotgarch0$results$par[2])-1)
  VaR95_pot_z0=quantile(ztm0,0.9)+gpdpotgarch0$results$par[1]/gpdpotgarch0$results$par[2]*((2400*0.05/240)^(-gpdpotgarch0$results$par[2])-1)
  
  VaR995_pot_garch0[i]=VaR995_pot_z0*predict(garchfitm0)[1,3]+predict(garchfitm0)[1,1]#Mcneil S.6
  VaR99_pot_garch0[i]=VaR99_pot_z0*predict(garchfitm0)[1,3]+predict(garchfitm0)[1,1]
  VaR95_pot_garch0[i]=VaR95_pot_z0*predict(garchfitm0)[1,3]+predict(garchfitm0)[1,1]
  #ES: Mcneil S.293
  ES995_pot_garch0[i]=predict(garchfitm0)[1,1]+predict(garchfitm0)[1,3]*VaR995_pot_z0*(1/(1-gpdpotgarch0$results$par[2])+(gpdpotgarch0$results$par[1]-gpdpotgarch0$results$par[2]*quantile(ztm0,0.9))/(VaR995_pot_z0-VaR995_pot_z0*gpdpotgarch0$results$par[2]))
  ES99_pot_garch0[i]=predict(garchfitm0)[1,1]+predict(garchfitm0)[1,3]*VaR99_pot_z0*(1/(1-gpdpotgarch0$results$par[2])+(gpdpotgarch0$results$par[1]-gpdpotgarch0$results$par[2]*quantile(ztm0,0.9))/(VaR99_pot_z0-VaR995_pot_z0*gpdpotgarch0$results$par[2]))
  ES95_pot_garch0[i]=predict(garchfitm0)[1,1]+predict(garchfitm0)[1,3]*VaR95_pot_z0*(1/(1-gpdpotgarch0$results$par[2])+(gpdpotgarch0$results$par[1]-gpdpotgarch0$results$par[2]*quantile(ztm0,0.9))/(VaR95_pot_z0-VaR995_pot_z0*gpdpotgarch0$results$par[2]))
}  
VaR995_pot_xts_garch0=xts(VaR995_pot_garch0,dax_log$date[2401:6826])
VaR99_pot_xts_garch0=xts(VaR99_pot_garch0,dax_log$date[2401:6826])
VaR95_pot_xts_garch0=xts(VaR95_pot_garch0,dax_log$date[2401:6826])

ES995_pot_xts_garch0=xts(ES995_pot_garch0,dax_log$date[2401:6826])
ES99_pot_xts_garch0=xts(ES99_pot_garch0,dax_log$date[2401:6826])
ES95_pot_xts_garch0=xts(ES95_pot_garch0,dax_log$date[2401:6826])

plot(dax_log_xts[2401:6826],main="VaR",ylim=c(0,12))  
lines(VaR995_pot_xts_garch0,col="red")   
lines(VaR99_pot_xts_garch0,col="blue")    
lines(VaR95_pot_xts_garch0,col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_pot_xts_garch0<dax_log_xts[2401:6826]) #  16 Ueberschreitungen
sum(VaR99_pot_xts_garch0<dax_log_xts[2401:6826]) #  35 Ueberschreitungen
sum(VaR95_pot_xts_garch0<dax_log_xts[2401:6826]) #  232 Ueberschreitungen

plot(dax_log_xts[2401:6826],main="ES",ylim=c(0,12))  
lines(ES995_pot_xts_garch0,col="red")   
lines(ES99_pot_xts_garch0,col="blue")    
lines(ES95_pot_xts_garch0,col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#nur die ersten 1000 Prognosen
plot(dax_log_xts[2401:3400],main="VaR",ylim=c(0,12))  
lines(VaR995_pot_xts_garch0[1:1000],col="red")   
lines(VaR99_pot_xts_garch0[1:1000],col="blue")    
lines(VaR95_pot_xts_garch0[1:1000],col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

plot(dax_log_xts[2401:3400],main="ES",ylim=c(0,13))  
lines(ES995_pot_xts_garch0[1:1000],col="red")   
lines(ES99_pot_xts_garch0[1:1000],col="blue")    
lines(ES95_pot_xts_garch0[1:1000],col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_pot_xts_garch0<dax_log_xts[2401:6826])
sum(V995)
V99=(VaR99_pot_xts_garch0<dax_log_xts[2401:6826])
sum(V99)
V95=(VaR95_pot_xts_garch0<dax_log_xts[2401:6826])
sum(V95)

uc_test(p=0.005,v=V995) #1.889475
uc_test(p=0.01,v=V99)  #2.108226
uc_test(p=0.05,v=V95)  #0.536462
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.1692612
1-pchisq(uc_test(p=0.01,v=V99),1)#0.146509
1-pchisq(uc_test(p=0.05,v=V95),1)#0.4639027

##Ind.Test 
ind_test(as.vector(V995))#3.929223
ind_test(as.vector(V99)) #1.152785
ind_test(as.vector(V95))  #0.06278441
1-pchisq(ind_test(as.vector(V995)),1)#0.04745386
1-pchisq(ind_test(as.vector(V99)),1)#0.282967
1-pchisq(ind_test(as.vector(V95)),1)#0.802148

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.0545112
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.1958306
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.7410974


####Acerbi Test 2
Z2=function(p,ES,L,v){
  s = matrix(ncol = 1,nrow = length(ES))
  for (i in 1:length(ES)){
    s[i]=L[i]*v[i]/(p*length(ES)*ES[i])
  }
  return(sum(s)-1)
}
Z2(p=0.005,ES=ES995_pot_xts_garch0,L=dax_log_xts[2401:6826],v=V995)#-0.1895934
Z2(p=0.01,ES=ES99_pot_xts_garch0,L=dax_log_xts[2401:6826],v=V99)#-0.1590294
Z2(p=0.05,ES=ES95_pot_xts_garch0,L=dax_log_xts[2401:6826],v=V95)#0.04939282

esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES995_pot_xts_garch0,alpha=0.005)#0.1878786
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES99_pot_xts_garch0,alpha=0.01)#0.3548115
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES95_pot_xts_garch0,alpha=0.05)#0.3068557

###HS
sp=numeric(0)
VaR95_hs=numeric(0)
VaR99_hs=numeric(0)
VaR995_hs=numeric(0)
ES95_hs=numeric(0)
ES99_hs=numeric(0)
ES995_hs=numeric(0)

for (i in 1:4426){
  sp=dax_log_xts[i:(2400+i)]
  VaR95_hs[i]=quantile(sp,0.95)
  VaR99_hs[i]=quantile(sp,0.99)
  VaR995_hs[i]=quantile(sp,0.995)
  ES95_hs[i]=mean(sp[sp>VaR95_hs[i]])
  ES99_hs[i]=mean(sp[sp>VaR99_hs[i]])
  ES995_hs[i]=mean(sp[sp>VaR995_hs[i]])
}

plot(dax_log_xts[2401:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_hs,dax_log$date[2401:6826]),col="red")   
lines(xts(VaR99_hs,dax_log$date[2401:6826]),col="blue")    
lines(xts(VaR95_hs,dax_log$date[2401:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_hs<dax_log_xts[2401:6826]) #  23 Ueberschreitungen
sum(VaR99_hs<dax_log_xts[2401:6826]) #  51 Ueberschreitungen
sum(VaR95_hs<dax_log_xts[2401:6826]) #  218 Ueberschreitungen

plot(dax_log_xts[2401:6826],main="ES",ylim=c(0,10))  
lines(xts(ES995_hs,dax_log$date[2401:6826]),col="red")   
lines(xts(ES99_hs,dax_log$date[2401:6826]),col="blue")    
lines(xts(ES95_hs,dax_log$date[2401:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_hs<dax_log_xts[2401:6826])
sum(V995)
V99=(VaR99_hs<dax_log_xts[2401:6826])
sum(V99)
V95=(VaR95_hs<dax_log_xts[2401:6826])
sum(V95)

uc_test(p=0.005,v=V995) #0.03393473
uc_test(p=0.01,v=V99)  #0.9882912
uc_test(p=0.05,v=V95)  #0.05204494
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.8538457
1-pchisq(uc_test(p=0.01,v=V99),1)#0.3201604
1-pchisq(uc_test(p=0.05,v=V95),1)#0.8195424

##Ind.Test 
ind_test(as.vector(V995))#2.555994
ind_test(as.vector(V99)) #8.992341
ind_test(as.vector(V95))  #29.40033
1-pchisq(ind_test(as.vector(V995)),1)#0.1098767
1-pchisq(ind_test(as.vector(V99)),1)#0.002711134
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.2739076
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.006803512
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00


####Acerbi Test 2
Z2=function(p,ES,L,v){
  s = matrix(ncol = 1,nrow = length(ES))
  for (i in 1:length(ES)){
    s[i]=L[i]*v[i]/(p*length(ES)*ES[i])
  }
  return(sum(s)-1)
}
Z2(p=0.005,ES=ES995_hs,L=dax_log_xts[2401:6826],v=V995)#0.0420184
Z2(p=0.01,ES=ES99_hs,L=dax_log_xts[2401:6826],v=V99)#0.1416993
Z2(p=0.05,ES=ES95_hs,L=dax_log_xts[2401:6826],v=V95)#0.00695256

esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES995_hs,alpha=0.005)#0.4350271
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES99_hs,alpha=0.01)#0.3316457
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES95_hs,alpha=0.05)#0.2723359

##NV
sp=numeric(0)
VaR95_nv=numeric(0)
VaR99_nv=numeric(0)
VaR995_nv=numeric(0)
ES95_nv=numeric(0)
ES99_nv=numeric(0)
ES995_nv=numeric(0)

for (i in 1:4426){
  sp=dax_log_xts[i:(2400+i)]
  VaR95_nv[i]=mean(sp)+sd(sp)*qnorm(0.95)
  VaR99_nv[i]=mean(sp)+sd(sp)*qnorm(0.99)
  VaR995_nv[i]=mean(sp)+sd(sp)*qnorm(0.995)
  ES95_nv[i]=mean(sp)+sd(sp)*dnorm(qnorm(0.95))/0.05
  ES99_nv[i]=mean(sp)+sd(sp)*dnorm(qnorm(0.99))/0.01
  ES995_nv[i]=mean(sp)+sd(sp)*dnorm(qnorm(0.995))/0.005
}

plot(dax_log_xts[2401:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_nv,dax_log$date[2401:6826]),col="red")   
lines(xts(VaR99_nv,dax_log$date[2401:6826]),col="blue")    
lines(xts(VaR95_nv,dax_log$date[2401:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_nv<dax_log_xts[2401:6826]) #  70 Ueberschreitungen
sum(VaR99_nv<dax_log_xts[2401:6826]) #  87 Ueberschreitungen
sum(VaR95_nv<dax_log_xts[2401:6826]) #  218 Ueberschreitungen

plot(dax_log_xts[2401:6826],main="ES",ylim=c(0,10))  
lines(xts(ES995_nv,dax_log$date[2401:6826]),col="red")   
lines(xts(ES99_nv,dax_log$date[2401:6826]),col="blue")    
lines(xts(ES95_nv,dax_log$date[2401:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_nv<dax_log_xts[2401:6826])
sum(V995)
V99=(VaR99_nv<dax_log_xts[2401:6826])
sum(V99)
V95=(VaR95_nv<dax_log_xts[2401:6826])
sum(V95)

uc_test(p=0.005,v=V995) #66
uc_test(p=0.01,v=V99)  #32.53211
uc_test(p=0.05,v=V95)  #0.02524534
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0
1-pchisq(uc_test(p=0.01,v=V99),1)#0
1-pchisq(uc_test(p=0.05,v=V95),1)#0.8737573

##Ind.Test 
ind_test(as.vector(V995))#11.20981
ind_test(as.vector(V99)) #20.4193
ind_test(as.vector(V95))  #28.9707
1-pchisq(ind_test(as.vector(V995)),1)#0.00
1-pchisq(ind_test(as.vector(V99)),1)#0.00
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0


####Acerbi Test 2
Z2=function(p,ES,L,v){
  s = matrix(ncol = 1,nrow = length(ES))
  for (i in 1:length(ES)){
    s[i]=L[i]*v[i]/(p*length(ES)*ES[i])
  }
  return(sum(s)-1)
}
Z2(p=0.005,ES=ES995_nv,L=dax_log_xts[2401:6826],v=V995)#2.728504
Z2(p=0.01,ES=ES99_nv,L=dax_log_xts[2401:6826],v=V99)#1.379741
Z2(p=0.05,ES=ES95_nv,L=dax_log_xts[2401:6826],v=V95)#0.1736356

esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES995_nv,alpha=0.005)#0
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES99_nv,alpha=0.01)#0
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES95_nv,alpha=0.05)#0


##t
sp=numeric(0)
VaR95_t=numeric(0)
VaR99_t=numeric(0)
VaR995_t=numeric(0)
ES95_t=numeric(0)
ES99_t=numeric(0)
ES995_t=numeric(0)

for (i in 1:4426){
  sp=dax_log_xts[i:(2400+i)]
  tfit=fitdistr(sp,"t",lower=c(-1, 0.001, 1))
  para=as.numeric(tfit$estimate)
  VaR95_t[i]=para[1]+para[2]*qt(0.95,df=para[3])
  VaR99_t[i]=para[1]+para[2]*qt(0.99,df=para[3])
  VaR995_t[i]=para[1]+para[2]*qt(0.995,df=para[3])
  ES95_t[i]=para[1]+para[2]*(dt(qt(0.95,df=para[3]),df=para[3])/0.05*((para[3]+(qt(0.95,df=para[3]))^2)/(para[3]-1)))
  ES99_t[i]=para[1]+para[2]*(dt(qt(0.99,df=para[3]),df=para[3])/0.01*((para[3]+(qt(0.99,df=para[3]))^2)/(para[3]-1)))
  ES995_t[i]=para[1]+para[2]*(dt(qt(0.995,df=para[3]),df=para[3])/0.005*((para[3]+(qt(0.995,df=para[3]))^2)/(para[3]-1)))
}

plot(dax_log_xts[2401:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_t,dax_log$date[2401:6826]),col="red")   
lines(xts(VaR99_t,dax_log$date[2401:6826]),col="blue")    
lines(xts(VaR95_t,dax_log$date[2401:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_t<dax_log_xts[2401:6826]) #  25 Ueberschreitungen
sum(VaR99_t<dax_log_xts[2401:6826]) #  65 Ueberschreitungen
sum(VaR95_t<dax_log_xts[2401:6826]) #  265 Ueberschreitungen

plot(dax_log_xts[2401:6826],main="ES",ylim=c(0,12))  
lines(xts(ES995_t,dax_log$date[2401:6826]),col="red")   
lines(xts(ES99_t,dax_log$date[2401:6826]),col="blue")    
lines(xts(ES95_t,dax_log$date[2401:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_t<dax_log_xts[2401:6826])
sum(V995)
V99=(VaR99_t<dax_log_xts[2401:6826])
sum(V99)
V95=(VaR95_t<dax_log_xts[2401:6826])
sum(V95)

uc_test(p=0.005,v=V995) #0.3589543
uc_test(p=0.01,v=V99)  #8.578095
uc_test(p=0.05,v=V95)  #8.567368
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.5490876
1-pchisq(uc_test(p=0.01,v=V99),1)#0
1-pchisq(uc_test(p=0.05,v=V95),1)#0

##Ind.Test 
ind_test(as.vector(V995))#7.170128
ind_test(as.vector(V99)) #5.666412
ind_test(as.vector(V95))  #28.44295
1-pchisq(ind_test(as.vector(V995)),1)#0.007412749
1-pchisq(ind_test(as.vector(V99)),1)#0.0172928
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.02317825
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.00
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00


####Acerbi Test 2
Z2=function(p,ES,L,v){
  s = matrix(ncol = 1,nrow = length(ES))
  for (i in 1:length(ES)){
    s[i]=L[i]*v[i]/(p*length(ES)*ES[i])
  }
  return(sum(s)-1)
}
Z2(p=0.005,ES=ES995_t,L=dax_log_xts[2401:6826],v=V995)#0.03
Z2(p=0.01,ES=ES99_t,L=dax_log_xts[2401:6826],v=V99)#0.31
Z2(p=0.05,ES=ES95_t,L=dax_log_xts[2401:6826],v=V95)#0.17

esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES995_t,alpha=0.005)#0.9461494
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES99_t,alpha=0.01)#0.7387873
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES95_t,alpha=0.05)#0.1002435

#skt
sp=numeric(0)
VaR95_skt=numeric(0)
VaR99_skt=numeric(0)
VaR995_skt=numeric(0)
ES95_skt=numeric(0)
ES99_skt=numeric(0)
ES995_skt=numeric(0)

for (i in 1:4426){
  sp=dax_log_xts[i:(2400+i)]
  sktfit=sstdFit(sp)
  para=as.numeric(tfit$estimate)
  VaR95_t[i]=para[1]+para[2]*qt(0.95,df=para[3])
  VaR99_t[i]=para[1]+para[2]*qt(0.99,df=para[3])
  VaR995_t[i]=para[1]+para[2]*qt(0.995,df=para[3])
  ES95_t[i]=para[1]+para[2]*(dt(qt(0.95,df=para[3]),df=para[3])/0.05*((para[3]+(qt(0.95,df=para[3]))^2)/(para[3]-1)))
  ES99_t[i]=para[1]+para[2]*(dt(qt(0.99,df=para[3]),df=para[3])/0.01*((para[3]+(qt(0.99,df=para[3]))^2)/(para[3]-1)))
  ES995_t[i]=para[1]+para[2]*(dt(qt(0.995,df=para[3]),df=para[3])/0.005*((para[3]+(qt(0.995,df=para[3]))^2)/(para[3]-1)))
}

}

plot(dax_log_xts[2401:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_skt,dax_log$date[2401:6826]),col="red")   
lines(xts(VaR99_skt,dax_log$date[2401:6826]),col="blue")    
lines(xts(VaR95_skt,dax_log$date[2401:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_skt<dax_log_xts[2401:6826]) #  70 Ueberschreitungen
sum(VaR99_skt<dax_log_xts[2401:6826]) #  87 Ueberschreitungen
sum(VaR95_skt<dax_log_xts[2401:6826]) #  218 Ueberschreitungen

plot(dax_log_xts[2401:6826],main="ES",ylim=c(0,10))  
lines(xts(ES995_skt,dax_log$date[2401:6826]),col="red")   
lines(xts(ES99_skt,dax_log$date[2401:6826]),col="blue")    
lines(xts(ES95_skt,dax_log$date[2401:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_skt<dax_log_xts[2401:6826])
sum(V995)
V99=(VaR99_skt<dax_log_xts[2401:6826])
sum(V99)
V95=(VaR95_skt<dax_log_xts[2401:6826])
sum(V95)

uc_test(p=0.005,v=V995) #66
uc_test(p=0.01,v=V99)  #32.53211
uc_test(p=0.05,v=V95)  #0.02524534
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0
1-pchisq(uc_test(p=0.01,v=V99),1)#0
1-pchisq(uc_test(p=0.05,v=V95),1)#0.8737573

##Ind.Test 
ind_test(as.vector(V995))#11.20981
ind_test(as.vector(V99)) #20.4193
ind_test(as.vector(V95))  #28.9707
1-pchisq(ind_test(as.vector(V995)),1)#0.00
1-pchisq(ind_test(as.vector(V99)),1)#0.00
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0


####Acerbi Test 2
Z2=function(p,ES,L,v){
  s = matrix(ncol = 1,nrow = length(ES))
  for (i in 1:length(ES)){
    s[i]=L[i]*v[i]/(p*length(ES)*ES[i])
  }
  return(sum(s)-1)
}
Z2(p=0.005,ES=ES995_skt,L=dax_log_xts[2401:6826],v=V995)#2.728504
Z2(p=0.01,ES=ES99_skt,L=dax_log_xts[2401:6826],v=V99)#1.379741
Z2(p=0.05,ES=ES95_skt,L=dax_log_xts[2401:6826],v=V95)#0.1736356

esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES995_skt,alpha=0.005)#0
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES99_skt,alpha=0.01)#0
esr_backtest_intercept(-dax_log_xts[2401:6826],e=-ES95_skt,alpha=0.05)#0


##GARCH


