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
library(fGarch) #garch
library("vrtest") #Variance Ratio Test

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
length(jahrmax_bsp) #27 Jahre aber 28 Beobachtungen
plot(jahrmax_bsp)

# Maximum-Likelihood von GEV (extRemes Package)
fit_bsp <- fevd(as.vector(jahrmax_bsp), method = "MLE", type="GEV")
# Diagnostik Plots
plot(fit_bsp)
fit_bsp$results$par     #Paramter. Location=3.55 (mu), Scale=1.39 (sigma), Shape=0.195 (xi) #xi >0, somit ist es Frechet-Verteilung
# oder mit evir Package:
fit_bsp_evir=evir::gev(dax_log_xts,block = 253)  #1 Jahre hat etw. 253 Beobachtungen
plot(fit_bsp_evir) #gut gefittet
fit_bsp_evir$par.ests   #mu=3.6211049, sigma=1.4509886,xi=0.16


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

period.max(dax_log_xts,seq(from=21,to=6804,by=21))  #Monatlich z.B.

###Goodness of Fit Shermann 1957
bmm_monat_st=period.max(dax_log_xts[1:6804],seq(from=21,to=6804,by=21))
bmm_quartal_st=period.max(dax_log_xts[1:6804],seq(from=63,to=6804,by=63))
bmm_semester_st=period.max(dax_log_xts[1:6804],seq(from=126,to=6804,by=126))
bmm_jahr_st=period.max(dax_log_xts[1:6804],seq(from=252,to=6804,by=252))

####Monat:
fit_monat_st <- fevd(as.vector(bmm_monat_st), method = "MLE", type="GEV") 
x_monat_st=sort(fit_monat_st$x)    ####von klein zu gross geordnet
z_monat_st=pgev(x_monat_st,xi=fit_monat_st$results$par[3],mu=fit_monat_st$results$par[1],sigma=fit_monat_st$results$par[2])

###Goodness of Fit Shermann 1957
monat_order=sort(as.numeric(bmm_monat_st))
length(monat_order) #324
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
omega #statistik = 0.3830586
mean_monat=(324/325)^(324+1)
var_monat=(2*exp(1)-5)/(exp(2)*324)
qnorm(0.95,mean=mean_monat,sd=sqrt(var_monat)) #Kritischer Wert= 0.3895246 > 0.3830586 H0 nicht abgelehnt
pnorm(omega,mean=mean_monat,sd=sqrt(var_monat)) #p= 0.88

#omega1 oder 2
omega1=1/sqrt(var_monat)*(omega-mean_monat)
omega2=sqrt(exp(2)/(2*exp(1)-1)*324)*(omega-exp(-1))
omega1 #1.17
omega2 #0.35

###Quartal
fit_quartal_st <- fevd(as.vector(bmm_quartal_st), method = "MLE", type="GEV") 
x_quartal_st=sort(fit_quartal_st$x)    ####von klein zu gross geordnet
z_quartal_st=pgev(x_quartal_st,xi=fit_quartal_st$results$par[3],mu=fit_quartal_st$results$par[1],sigma=fit_quartal_st$results$par[2])
quartal_order=sort(as.numeric(bmm_quartal_st))
length(quartal_order) #108
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
omega #statistik = 0.3844832
mean_quartal=(108/109)^(108+1)
var_quartal=(2*exp(1)-5)/(exp(2)*108)
qnorm(0.95,mean=mean_quartal,sd=sqrt(var_quartal)) #Kritischer Wert= 0.4046574> 0.3844832 H0 nicht abgelehnt
pnorm(omega,mean=mean_quartal,sd=sqrt(var_quartal))  #p=0.7829847

#Semester
fit_semester_st <- fevd(as.vector(bmm_semester_st), method = "MLE", type="GEV") 
x_semester_st=sort(fit_semester_st$x)    ####von klein zu gross geordnet
z_semester_st=pgev(x_semester_st,xi=fit_semester_st$results$par[3],mu=fit_semester_st$results$par[1],sigma=fit_semester_st$results$par[2])
semester_order=sort(as.numeric(bmm_semester_st))
length(semester_order) #54
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
omega #statistik = 0.3707199
mean_semester=(54/55)^(54+1)
var_semester=(2*exp(1)-5)/(exp(2)*54)
qnorm(0.95,mean=mean_semester,sd=sqrt(var_semester)) #Kritischer Wert= 0.4189171>0.3707199 H0 nicht ab
pnorm(omega,mean=mean_semester,sd=sqrt(var_semester)) #0.5744649

##Jahr
fit_jahr_st <- fevd(as.vector(bmm_jahr_st), method = "MLE", type="GEV") 
x_jahr_st=sort(fit_jahr_st$x)    ####von klein zu gross geordnet
z_jahr_st=pgev(x_jahr_st,xi=fit_jahr_st$results$par[3],mu=fit_jahr_st$results$par[1],sigma=fit_jahr_st$results$par[2])
jahr_order=sort(as.numeric(bmm_jahr_st))
length(jahr_order) #27
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
omega #statistik = 0.3802239
mean_jahr=(27/28)^(27+1)
var_jahr=(2*exp(1)-5)/(exp(2)*27)
qnorm(0.95,mean=mean_jahr,sd=sqrt(var_jahr)) #Kritischer Wert= 0.4383879 > 0.3670251 H0 nicht ab.
pnorm(omega,mean=mean_jahr,sd=sqrt(var_jahr)) #0.66

####Alle bestehen den Test. Deshalb wird n = 21 gewaehlt! (mehr Beobachtungen) 

# Moving Window (Groesse=1050) Zum Beispiel: das erste Fenster. Laenge=1050, damit durch 21 perfekt teilbar.
ts_bm_1=dax_log_xts[1:1050]   
monatmax1=period.max(ts_bm_1,seq(from=21,to=1050,by=21))   #die groesseten monatlichen Verluste
monatmax1
plot(monatmax1)
fit_monat1 <- fevd(as.vector(monatmax1), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_monat1_evir=gev(dax_log_xts,block = 21)
plot(fit_monat1_evir)         #passt
fit_monat1$results$par  #Parameter extrahieren 

##fit an die ganzen Daten
monatmaxg=period.max(dax_log_xts,seq(from=21,to=6826,by=21))   #die groesseten monatlichen Verluste
monatmaxg
plot(monatmaxg)
fit_monatg <- fevd(as.vector(monatmaxg), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_monatg_evir=gev(dax_log_xts,block = 21)
fit_monatg_evir$par.ests
fit_monatg$results$par  #Parameter extrahieren 


# VaRs und ESs berechnen. n=21 ist monatlich
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
ES95_bmm=numeric(0)
ES99_bmm=numeric(0)
ES995_bmm=numeric(0)
VaR20 = function(x){mu-sigma/tau*(1-(-log((1-x)^(21)))^(-tau))} #tau = xi = shape Parameter
fit=numeric(0)
#returnlevel = function(x,mu,sigma,xi){mu-sigma/xi*(1-(-log((1-x)^21))^(-xi))}
for (i in (1:5776)){         #es gibt (6826-1050) Vorhersagen
  monatmax=period.max(dax_log_xts[i:(i+1049)],seq(from=21,to=1050,by=21))    #die groesseten monatlichen Verluste
  fit <- fevd(as.vector(monatmax), method = "MLE", type="GEV")
  #fit = gev(dax_log_xts[i:(i+1049)],block=120)
  #VaR995_bmm[i]=returnlevel(0.005,mu=fit$par.ests[3],sigma=fit$par.ests[2],xi=fit$par.ests[1])#rlevel.gev(fit,2.212322)[2]
   # VaR99_bmm[i]=returnlevel(0.01,mu=fit$par.ests[3],sigma=fit$par.ests[2],xi=fit$par.ests[1])#rlevel.gev(fit,1.427308)[2]
  #  VaR95_bmm[i]=returnlevel(0.05,mu=fit$par.ests[3],sigma=fit$par.ests[2],xi=fit$par.ests[1])#rlevel.gev(fit,1.002127)[2]
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 10.008750)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 21
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 5.255630)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.516442)[1]
   mu=as.numeric(fit$results$par[1])
    sigma=as.numeric(fit$results$par[2])
    tau=as.numeric(fit$results$par[3])
    ES995_bmm[i]=integrate(VaR20,lower=0,upper=0.005)$value*200
    ES99_bmm[i]=integrate(VaR20,lower=0,upper=0.01)$value*100
    ES95_bmm[i]=integrate(VaR20,lower=0,upper=0.05)$value*20
    }

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[1051:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[1051:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[1051:6826]) 
ES995_bmm_xts=xts(ES995_bmm,dax_log$date[1051:6826])
ES99_bmm_xts=xts(ES99_bmm,dax_log$date[1051:6826])
ES95_bmm_xts=xts(ES95_bmm,dax_log$date[1051:6826])

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))


plot(dax_log_xts[1051:6826],main="ES",,ylim=c(0,10))  
lines(ES995_bmm_xts,col="red")   
lines(ES99_bmm_xts,col="blue")    
lines(ES95_bmm_xts,col="green")  
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_bmm_xts<dax_log_xts[1051:6826]) # 63 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[1051:6826]) #  121 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[1051:6826]) #  505 Ueberschreitungen

#U.C Test (Code aus dem Buch von Danielsson)
uc_test = function(p,v){
  return(2*log(((1-sum(v)/length(v))/(1-p))^(length(v)-sum(v))*((sum(v)/length(v))/p)^(sum(v))))
}

V995=(VaR995_bmm_xts<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_bmm_xts<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_bmm_xts<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #30.24113
uc_test(p=0.01,v=V99)  #53.17955
uc_test(p=0.05,v=V95)  #140.6448 
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1) #0
1-pchisq(uc_test(p=0.01,v=V99),1) #0
1-pchisq(uc_test(p=0.05,v=V95),1) #0

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

ind_test(as.vector(V995))#4.391147
ind_test(as.vector(V99)) #16.61035
ind_test(as.vector(V95))  #44.50534
1-pchisq(ind_test(as.vector(V995)),1)#0.03612599
1-pchisq(ind_test(as.vector(V99)),1)#0.00
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.00
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.00
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00

##ES Test (Mcneil 2000)
#ESTest(0.005,-dax_log_xts[1051:6826],ES=-ES995_bmm_xts,VaR=-VaR995_bmm_xts)
#ESTest(0.01,-dax_log_xts[1051:6826],ES=-ES99_bmm_xts,VaR=-VaR99_bmm_xts)
#ESTest(0.05,actual=-dax_log_xts[1051:6826],ES=-ES95_bmm_xts,VaR=-VaR95_bmm_xts)

####Acerbi Test 2
Z2=function(p,ES,L,v){
  s = matrix(ncol = 1,nrow = length(ES))
  for (i in 1:length(ES)){
  s[i]=L[i]*v[i]/(p*length(ES)*ES[i])
  }
  return(sum(s)-1)
}
Z2(p=0.005,ES=ES995_bmm_xts,L=dax_log_xts[1051:6826],v=V995)#1.220596
Z2(p=0.01,ES=ES99_bmm_xts,L=dax_log_xts[1051:6826],v=V99)#1.127528
Z2(p=0.05,ES=ES95_bmm_xts,L=dax_log_xts[1051:6826],v=V95)#0.8269118

#ESBACK
#esback::esr_backtest(-dax_log_xts[1051:6826],e=-ES995_bmm_xts,alpha=0.005)
#esback::esr_backtest(-dax_log_xts[1051:6826],e=-ES99_bmm_xts,alpha=0.01)
#esback::esr_backtest(-dax_log_xts[1051:6826],e=-ES95_bmm_xts,alpha=0.05)

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_bmm_xts,alpha=0.005)#0.00
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_bmm_xts,alpha=0.01)#0.00
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_bmm_xts,alpha=0.05)#0.00

#esback::cc_backtest(-dax_log_xts[1051:6826],q=-ES995_bmm_xts,e=-VaR995_bmm_xts,alpha=0.005)
#esback::cc_backtest(-dax_log_xts[1051:6826],q=-ES99_bmm_xts,e=-VaR99_bmm_xts,alpha=0.01)
#esback::cc_backtest(-dax_log_xts[1051:6826],q=-ES99_bmm_xts,e=-VaR95_bmm_xts,alpha=0.05)

#esback::er_backtest(-dax_log_xts[1051:6826],q=-ES995_bmm_xts,e=-VaR995_bmm_xts)
#esback::er_backtest(-dax_log_xts[1051:6826],q=-ES99_bmm_xts,e=-VaR99_bmm_xts)
#esback::er_backtest(-dax_log_xts[1051:6826],q=-ES95_bmm_xts,e=-VaR95_bmm_xts)



# Mit Theta (Extremaler Index) Embrechts 1998 Chap8 P.S19

#n=63,n=108 N_u=15,20,25,30,40,50,100,200

#n=63  #Mcneil 1998 S.13,14  (m, n muessen gross genug sein)
n60max=period.max(dax_log_xts[1:6804],seq(from=63,to=6804,by=63))
fit_60 <- fevd(as.vector(n60max), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
#plot(fit_60)
#return.level(fit_60, conf = 0.05, return.period= c(2,5,10,20,50))
N15=quantile(dax_log_xts,(6804-15)/6804) #N_u=15
N20=quantile(dax_log_xts,(6804-20)/6804) #N_u=20
N25=quantile(dax_log_xts,(6804-25)/6804) #N_u=25
N30=quantile(dax_log_xts,(6804-30)/6804) #N_u=30
N40=quantile(dax_log_xts,(6804-40)/6804) #N_u=40
N50=quantile(dax_log_xts,(6804-50)/6804) #N_u=50
N100=quantile(dax_log_xts,(6804-100)/6804) #N_u=100
N200=quantile(dax_log_xts,(6804-200)/6804) #N_u=200
K15=sum(n60max>N15);K20=sum(n60max>N20);K25=sum(n60max>N25);K30=sum(n60max>N30);K40=sum(n60max>N40);
K50=sum(n60max>N50);K100=sum(n60max>N100);K200=sum(n60max>N200)
K15;K20;K25;K30;K40;K50;K100;K200
theta=function(n,m,Ku,Nu){(log(1-Ku/m)/log(1-Nu/6804))/n}
#Theta berechnen
theta15=theta(n=63,m=108,Ku=K15,Nu=15)
theta20=theta(n=63,m=108,Ku=K20,Nu=20)
theta25=theta(n=63,m=108,Ku=K25,Nu=25)
theta30=theta(n=63,m=108,Ku=K30,Nu=30)
theta40=theta(n=63,m=108,Ku=K40,Nu=40)
theta50=theta(n=63,m=108,Ku=K50,Nu=50)
theta100=theta(n=63,m=108,Ku=K100,Nu=100)
theta200=theta(n=63,m=108,Ku=K200,Nu=200)
theta_dach=mean(c(theta15,theta20,theta25,theta30,theta40,theta50,theta100,theta200))#Mcneil 1998 S.13,14
theta_dach #0.52358 

#n=126
n120max=period.max(dax_log_xts[1:6804],seq(from=126,to=6804,by=126))
fit_120 <- fevd(as.vector(n120max), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
#plot(fit_60)
#return.level(fit_60, conf = 0.05, return.period= c(2,5,10,20,50))
N15=quantile(dax_log_xts,(6804-15)/6804) #N_u=15
N20=quantile(dax_log_xts,(6804-20)/6804) #N_u=20
N25=quantile(dax_log_xts,(6804-25)/6804) #N_u=25
N30=quantile(dax_log_xts,(6804-30)/6804) #N_u=30
N40=quantile(dax_log_xts,(6804-40)/6804) #N_u=40
N50=quantile(dax_log_xts,(6804-50)/6804) #N_u=50
N100=quantile(dax_log_xts,(6804-100)/6804) #N_u=100
N200=quantile(dax_log_xts,(6804-200)/6804) #N_u=200
K15=sum(n120max>N15);K20=sum(n120max>N20);K25=sum(n120max>N25);K30=sum(n120max>N30);K40=sum(n120max>N40);
K50=sum(n120max>N50);K100=sum(n120max>N100);K200=sum(n120max>N200)
K15;K20;K25;K30;K40;K50;K100;K200
theta=function(n,m,Ku,Nu){(log(1-Ku/m)/log(1-Nu/6804))/n}
#Theta berechnen
theta15=theta(n=126,m=54,Ku=K15,Nu=15)
theta20=theta(n=126,m=54,Ku=K20,Nu=20)
theta25=theta(n=126,m=54,Ku=K25,Nu=25)
theta30=theta(n=126,m=54,Ku=K30,Nu=30)
theta40=theta(n=126,m=54,Ku=K40,Nu=40)
theta50=theta(n=126,m=54,Ku=K50,Nu=50)
theta100=theta(n=126,m=54,Ku=K100,Nu=100)
theta200=theta(n=126,m=54,Ku=K200,Nu=200)
theta_dach=mean(c(theta15,theta20,theta25,theta30,theta40,theta50,theta100,theta200))#Mcneil 1998 S.13,14
theta_dach #0.44889  



#Theta mir "evir" berechnen. n=126 ...

#VaR
#n=126, Moving Window = 1050 Theta= 0.52. (S.13 Mcneil1998)
VaR95_bmme=numeric(0)
VaR99_bmme=numeric(0)
VaR995_bmme=numeric(0)
ES95_bmme=numeric(0)
ES99_bmme=numeric(0)
ES995_bmme=numeric(0)
n20max=numeric(0)
fit=numeric(0)
VaR20 = function(x){mu-sigma/tau*(1-(-log((1-x)^(21*0.52)))^(-tau))} #tau = xi = shape Parameter

for (i in (1:5776)){         #es gibt (6826-1050) Vorhersagen
  n20max=period.max(dax_log_xts[i:(i+1049)],seq(from=21,to=1050,by=21))    #der groessete quartalliche Verlust
  fit <- fevd(as.vector(n20max), method = "MLE", type="GEV")
  VaR995_bmme[i]=return.level(fit, conf = 0.05, return.period= 18.773754)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^(n*theta)=(1-1/k). Hier n = 21, Theta=0.52
  VaR99_bmme[i]=return.level(fit, conf = 0.05, return.period= 9.620789)[1]  
  VaR95_bmme[i]=return.level(fit, conf = 0.05, return.period= 2.3317576)[1] 
  mu=as.numeric(fit$results$par[1])
  sigma=as.numeric(fit$results$par[2])
  tau=as.numeric(fit$results$par[3])
  ES995_bmme[i]=integrate(VaR20,lower=0,upper=0.005)$value*200
  ES99_bmme[i]=integrate(VaR20,lower=0,upper=0.01)$value*100
  ES95_bmme[i]=integrate(VaR20,lower=0,upper=0.05)$value*20
}

VaR995_bmme_xts=xts(VaR995_bmme,dax_log$date[1051:6826])
VaR99_bmme_xts=xts(VaR99_bmme,dax_log$date[1051:6826])
VaR95_bmme_xts=xts(VaR95_bmme,dax_log$date[1051:6826]) 
ES995_bmme_xts=xts(ES995_bmme,dax_log$date[1051:6826])
ES99_bmme_xts=xts(ES99_bmme,dax_log$date[1051:6826])
ES95_bmme_xts=xts(ES95_bmme,dax_log$date[1051:6826])

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_bmme_xts,col="red")   
lines(VaR99_bmme_xts,col="blue")    
lines(VaR95_bmme_xts,col="green")  
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))


plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
lines(ES995_bmme_xts,col="red")   
lines(ES99_bmme_xts,col="blue")    
lines(ES95_bmme_xts,col="green")  
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_bmme_xts<dax_log_xts[1051:6826]) #  38 Ueberschreitungen
sum(VaR99_bmme_xts<dax_log_xts[1051:6826]) #  63 Ueberschreitungen
sum(VaR95_bmme_xts<dax_log_xts[1051:6826]) #  289 Ueberschreitungen

#U.C Test
V995=(VaR995_bmme_xts<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_bmme_xts<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_bmme_xts<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #2.63168
uc_test(p=0.01,v=V99)  #0.4664204
uc_test(p=0.05,v=V95)  #0.000145762
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.1047508
1-pchisq(uc_test(p=0.01,v=V99),1)#0.4946386
1-pchisq(uc_test(p=0.05,v=V95),1)#0.9903672

##Ind.Test 
ind_test(as.vector(V995))#1.302425
ind_test(as.vector(V99)) #4.391147
ind_test(as.vector(V95))  #32.93956
1-pchisq(ind_test(as.vector(V995)),1)#0.2537708
1-pchisq(ind_test(as.vector(V99)),1)#0.04
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.1398685
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.08814395
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00

####ES TEST

Z2(p=0.005,ES=ES995_bmme_xts,L=dax_log_xts[1051:6826],v=V995)#0.3274358
Z2(p=0.01,ES=ES99_bmme_xts,L=dax_log_xts[1051:6826],v=V99)#0.1216147
Z2(p=0.05,ES=ES95_bmme_xts,L=dax_log_xts[1051:6826],v=V95)#0.03152348

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_bmme_xts,alpha=0.005)#0.1163567
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_bmme_xts,alpha=0.01)#0.1354256
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_bmme_xts,alpha=0.05)#0.1592234


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


dax_log_order=sort(-dax_log$logreturn,decreasing = TRUE) 
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

#Moving Windows mit Laenge 1050
VaR95_pot=numeric(0)
VaR99_pot=numeric(0)
VaR995_pot=numeric(0)
ES95_pot=numeric(0)
ES99_pot=numeric(0)
ES995_pot=numeric(0)
for (i in (1:5776)){         #es gibt (6826-1050) Vorhersagen  
  gpdpot=gpd(dax_log_xts[i:(1049+i)],threshold=quantile(dax_log_xts[i:(1049+i)],0.9),method = "ml")
  #Mcneil,2000 S.288 the choice of k in Moving Window (10%)
  risk=riskmeasures(gpdpot,c(0.95,0.99,0.995))
  VaR995_pot[i]=risk[6]
  VaR99_pot[i]=risk[5]
  VaR95_pot[i]=risk[4]
  ES995_pot[i]=risk[9]
  ES99_pot[i]=risk[8]
  ES95_pot[i]=risk[7]
}  

VaR995_pot_xts=xts(VaR995_pot,dax_log$date[1051:6826])
VaR99_pot_xts=xts(VaR99_pot,dax_log$date[1051:6826])
VaR95_pot_xts=xts(VaR95_pot,dax_log$date[1051:6826]) 

ES995_pot_xts=xts(ES995_pot,dax_log$date[1051:6826])
ES99_pot_xts=xts(ES99_pot,dax_log$date[1051:6826])
ES95_pot_xts=xts(ES95_pot,dax_log$date[1051:6826])


plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_pot_xts,col="red")   
lines(VaR99_pot_xts,col="blue")    
lines(VaR95_pot_xts,col="green") 
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10))  
lines(ES995_pot_xts,col="red")   
lines(ES99_pot_xts,col="blue")    
lines(ES95_pot_xts,col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_pot_xts<dax_log_xts[1051:6826]) #  42 Ueberschreitungen
sum(VaR99_pot_xts<dax_log_xts[1051:6826]) #  73 Ueberschreitungen
sum(VaR95_pot_xts<dax_log_xts[1051:6826]) #  325 Ueberschreitungen

#U.C Test
V995=(VaR995_pot_xts<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_pot_xts<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_pot_xts<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #5.24968
uc_test(p=0.01,v=V99)  #3.748443
uc_test(p=0.05,v=V95)  #4.598348
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.02195081
1-pchisq(uc_test(p=0.01,v=V99),1)#0.05285672
1-pchisq(uc_test(p=0.05,v=V95),1)#0.03200277

##Ind.Test 
ind_test(as.vector(V995))#1.006171
ind_test(as.vector(V99)) #9.216849
ind_test(as.vector(V95))  #28.30808
1-pchisq(ind_test(as.vector(V995)),1)#0.3158218
1-pchisq(ind_test(as.vector(V99)),1)#0.002397979
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.04380858
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.001529758
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00


####ES TEST

Z2(p=0.005,ES=ES995_pot_xts,L=dax_log_xts[1051:6826],v=V995)#0.5674553
Z2(p=0.01,ES=ES99_pot_xts,L=dax_log_xts[1051:6826],v=V99)#0.3507733
Z2(p=0.05,ES=ES95_pot_xts,L=dax_log_xts[1051:6826],v=V95)#0.1842373

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_pot_xts,alpha=0.005)#0.008116015
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_pot_xts,alpha=0.01)#0.009651432
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_pot_xts,alpha=0.05)#0.001053563


##POT mit GARCH-Filter

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

#Moving Windows mit Laenge 1050 (Fuer jedes MovingWindow wird GARCH erneut geschaetzt)
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

for (i in (1:5776)){         #es gibt (6826-1050) Vorhersagen. (laeuft 45 Minuten..)
  garchfitm=garchFit(formula=~arma(1,0)+garch(1,1),data=dax_log_xts[i:(1049+i)],cond.dist ="QMLE")
  ztm=(dax_log_xts[i:(1049+i)]-garchfitm@fitted)/garchfitm@sigma.t
  gpdpotgarch=fevd(as.vector(ztm), method = "MLE", type="GP", threshold=quantile(ztm,(1050-105)/1050))
  VaR995_pot_z=quantile(ztm,0.9)+gpdpotgarch$results$par[1]/gpdpotgarch$results$par[2]*((1050*0.005/105)^(-gpdpotgarch$results$par[2])-1)
  VaR99_pot_z=quantile(ztm,0.9)+gpdpotgarch$results$par[1]/gpdpotgarch$results$par[2]*((1050*0.01/105)^(-gpdpotgarch$results$par[2])-1)
  VaR95_pot_z=quantile(ztm,0.9)+gpdpotgarch$results$par[1]/gpdpotgarch$results$par[2]*((1050*0.05/105)^(-gpdpotgarch$results$par[2])-1)
  
  VaR995_pot_garch[i]=VaR995_pot_z*predict(garchfitm)[1,3]+predict(garchfitm)[1,1]#Mcneil S.6
  VaR99_pot_garch[i]=VaR99_pot_z*predict(garchfitm)[1,3]+predict(garchfitm)[1,1]
  VaR95_pot_garch[i]=VaR95_pot_z*predict(garchfitm)[1,3]+predict(garchfitm)[1,1]
  #ES: Mcneil S.293
  ES995_pot_garch[i]=predict(garchfitm)[1,1]+predict(garchfitm)[1,3]*VaR995_pot_z*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(ztm,0.9))/(VaR995_pot_z-VaR995_pot_z*gpdpotgarch$results$par[2]))
  ES99_pot_garch[i]=predict(garchfitm)[1,1]+predict(garchfitm)[1,3]*VaR99_pot_z*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(ztm,0.9))/(VaR99_pot_z-VaR995_pot_z*gpdpotgarch$results$par[2]))
  ES95_pot_garch[i]=predict(garchfitm)[1,1]+predict(garchfitm)[1,3]*VaR95_pot_z*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(ztm,0.9))/(VaR95_pot_z-VaR995_pot_z*gpdpotgarch$results$par[2]))
}  
VaR995_pot_xts_garch=xts(VaR995_pot_garch,dax_log$date[1051:6826])
VaR99_pot_xts_garch=xts(VaR99_pot_garch,dax_log$date[1051:6826])
VaR95_pot_xts_garch=xts(VaR95_pot_garch,dax_log$date[1051:6826])

ES995_pot_xts_garch=xts(ES995_pot_garch,dax_log$date[1051:6826])
ES99_pot_xts_garch=xts(ES99_pot_garch,dax_log$date[1051:6826])
ES95_pot_xts_garch=xts(ES95_pot_garch,dax_log$date[1051:6826])

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_pot_xts_garch,col="red")   
lines(VaR99_pot_xts_garch,col="blue")    
lines(VaR95_pot_xts_garch,col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_pot_xts_garch<dax_log_xts[1051:6826]) #  26 Ueberschreitungen
sum(VaR99_pot_xts_garch<dax_log_xts[1051:6826]) #  52 Ueberschreitungen
sum(VaR95_pot_xts_garch<dax_log_xts[1051:6826]) #  308 Ueberschreitungen

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10))  
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


for (i in (1:5776)){         #es gibt (6826-1050) Vorhersagen. (laeuft 45 Minuten..)
  garchfitm0=garchFit(formula=~arma(0,0)+garch(1,1),data=dax_log_xts[i:(1049+i)],cond.dist ="QMLE")
  ztm0=(dax_log_xts[i:(1049+i)]-garchfitm0@fitted)/garchfitm0@sigma.t
  gpdpotgarch0=fevd(as.vector(ztm0), method = "MLE", type="GP", threshold=quantile(ztm0,(1050-105)/1050))
  VaR995_pot_z0=quantile(ztm0,0.9)+gpdpotgarch0$results$par[1]/gpdpotgarch0$results$par[2]*((1050*0.005/105)^(-gpdpotgarch0$results$par[2])-1)
  VaR99_pot_z0=quantile(ztm0,0.9)+gpdpotgarch0$results$par[1]/gpdpotgarch0$results$par[2]*((1050*0.01/105)^(-gpdpotgarch0$results$par[2])-1)
  VaR95_pot_z0=quantile(ztm0,0.9)+gpdpotgarch0$results$par[1]/gpdpotgarch0$results$par[2]*((1050*0.05/105)^(-gpdpotgarch0$results$par[2])-1)
  
  VaR995_pot_garch0[i]=VaR995_pot_z0*predict(garchfitm0)[1,3]+predict(garchfitm0)[1,1]#Mcneil S.6
  VaR99_pot_garch0[i]=VaR99_pot_z0*predict(garchfitm0)[1,3]+predict(garchfitm0)[1,1]
  VaR95_pot_garch0[i]=VaR95_pot_z0*predict(garchfitm0)[1,3]+predict(garchfitm0)[1,1]
  #ES: Mcneil S.293
  ES995_pot_garch0[i]=predict(garchfitm0)[1,1]+predict(garchfitm0)[1,3]*VaR995_pot_z0*(1/(1-gpdpotgarch0$results$par[2])+(gpdpotgarch0$results$par[1]-gpdpotgarch0$results$par[2]*quantile(ztm0,0.9))/(VaR995_pot_z0-VaR995_pot_z0*gpdpotgarch0$results$par[2]))
  ES99_pot_garch0[i]=predict(garchfitm0)[1,1]+predict(garchfitm0)[1,3]*VaR99_pot_z0*(1/(1-gpdpotgarch0$results$par[2])+(gpdpotgarch0$results$par[1]-gpdpotgarch0$results$par[2]*quantile(ztm0,0.9))/(VaR99_pot_z0-VaR995_pot_z0*gpdpotgarch0$results$par[2]))
  ES95_pot_garch0[i]=predict(garchfitm0)[1,1]+predict(garchfitm0)[1,3]*VaR95_pot_z0*(1/(1-gpdpotgarch0$results$par[2])+(gpdpotgarch0$results$par[1]-gpdpotgarch0$results$par[2]*quantile(ztm0,0.9))/(VaR95_pot_z0-VaR995_pot_z0*gpdpotgarch0$results$par[2]))
}  
VaR995_pot_xts_garch0=xts(VaR995_pot_garch0,dax_log$date[1051:6826])
VaR99_pot_xts_garch0=xts(VaR99_pot_garch0,dax_log$date[1051:6826])
VaR95_pot_xts_garch0=xts(VaR95_pot_garch0,dax_log$date[1051:6826])

ES995_pot_xts_garch0=xts(ES995_pot_garch0,dax_log$date[1051:6826])
ES99_pot_xts_garch0=xts(ES99_pot_garch0,dax_log$date[1051:6826])
ES95_pot_xts_garch0=xts(ES95_pot_garch0,dax_log$date[1051:6826])

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,12))  
lines(VaR995_pot_xts_garch0,col="red")   
lines(VaR99_pot_xts_garch0,col="blue")    
lines(VaR95_pot_xts_garch0,col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_pot_xts_garch0<dax_log_xts[1051:6826]) #  26 Ueberschreitungen
sum(VaR99_pot_xts_garch0<dax_log_xts[1051:6826]) #  51 Ueberschreitungen
sum(VaR95_pot_xts_garch0<dax_log_xts[1051:6826]) #  311 Ueberschreitungen

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,13))  
lines(ES995_pot_xts_garch0,col="red")   
lines(ES99_pot_xts_garch0,col="blue")    
lines(ES95_pot_xts_garch0,col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

##mit Punkten
plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10),type="p",pch=1)
cc=dax_log_xts[1051:6826]
points(cc[cc>VaR95_pot_garch0],type="p",col="blue",pch=16)
points(cc[cc>ES95_pot_garch0],type="p",col="red",pch=16)

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10),type="p",pch=1)
cc=dax_log_xts[1051:6826]
points(cc[cc>VaR99_pot_garch0],type="p",col="blue",pch=16)
points(cc[cc>ES99_pot_garch0],type="p",col="red",pch=16)

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10),type="p",pch=1)
cc=dax_log_xts[1051:6826]
points(cc[cc>VaR995_pot_garch0],type="p",col="blue",pch=16)
points(cc[cc>ES995_pot_garch0],type="p",col="red",pch=16)

#nur die ersten 1000 Prognosen
plot(dax_log_xts[1051:2050],main="VaR",ylim=c(0,12))  
lines(VaR995_pot_xts_garch0[1:1000],col="red")   
lines(VaR99_pot_xts_garch0[1:1000],col="blue")    
lines(VaR95_pot_xts_garch0[1:1000],col="green")   
legend("topleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

plot(dax_log_xts[1051:2050],main="ES",ylim=c(0,13))  
lines(ES995_pot_xts_garch0[1:1000],col="red")   
lines(ES99_pot_xts_garch0[1:1000],col="blue")    
lines(ES95_pot_xts_garch0[1:1000],col="green") 
legend("topleft",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_pot_xts_garch0<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_pot_xts_garch0<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_pot_xts_garch0<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #0.2986986
uc_test(p=0.01,v=V99)  #0.8319605
uc_test(p=0.05,v=V95)  #1.754328
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.5846994
1-pchisq(uc_test(p=0.01,v=V99),1)#0.3617062
1-pchisq(uc_test(p=0.05,v=V95),1)#0.1853336

##Ind.Test 
ind_test(as.vector(V995))#2.585386
ind_test(as.vector(V99)) #0.508112
ind_test(as.vector(V95))  #0.5305087
1-pchisq(ind_test(as.vector(V995)),1)#0.1078541
1-pchisq(ind_test(as.vector(V99)),1)#0.4759573
1-pchisq(ind_test(as.vector(V95)),1)#0.4663931

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.2364443
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.51169
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.3190465


###ES TEST

Z2(p=0.005,ES=ES995_pot_xts_garch0,L=dax_log_xts[1051:6826],v=V995)#-0.02996437
Z2(p=0.01,ES=ES99_pot_xts_garch0,L=dax_log_xts[1051:6826],v=V99)#-0.07229775
Z2(p=0.05,ES=ES95_pot_xts_garch0,L=dax_log_xts[1051:6826],v=V95)#0.07756284

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_pot_xts_garch0,alpha=0.005)#0.1245075
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_pot_xts_garch0,alpha=0.01)#0.2770953
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_pot_xts_garch0,alpha=0.05)#0.1570921

###HS
sp=numeric(0)
VaR95_hs=numeric(0)
VaR99_hs=numeric(0)
VaR995_hs=numeric(0)
ES95_hs=numeric(0)
ES99_hs=numeric(0)
ES995_hs=numeric(0)

for (i in 1:5776){
  sp=dax_log_xts[i:(1049+i)]
  VaR95_hs[i]=quantile(sp,0.95)
  VaR99_hs[i]=quantile(sp,0.99)
  VaR995_hs[i]=quantile(sp,0.995)
  ES95_hs[i]=mean(sp[sp>VaR95_hs[i]])
  ES99_hs[i]=mean(sp[sp>VaR99_hs[i]])
  ES995_hs[i]=mean(sp[sp>VaR995_hs[i]])
}

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_hs,dax_log$date[1051:6826]),col="red")   
lines(xts(VaR99_hs,dax_log$date[1051:6826]),col="blue")    
lines(xts(VaR95_hs,dax_log$date[1051:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_hs<dax_log_xts[1051:6826]) #  41 Ueberschreitungen
sum(VaR99_hs<dax_log_xts[1051:6826]) #  80 Ueberschreitungen
sum(VaR95_hs<dax_log_xts[1051:6826]) #  323 Ueberschreitungen

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10))  
lines(xts(ES995_hs,dax_log$date[1051:6826]),col="red")   
lines(xts(ES99_hs,dax_log$date[1051:6826]),col="blue")    
lines(xts(ES95_hs,dax_log$date[1051:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_hs<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_hs<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_hs<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #4.520243
uc_test(p=0.01,v=V99)  #7.723433
uc_test(p=0.05,v=V95)  #4.112576
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.03349607
1-pchisq(uc_test(p=0.01,v=V99),1)#0.005450865
1-pchisq(uc_test(p=0.05,v=V95),1)#0.04256548

##Ind.Test 
ind_test(as.vector(V995))#1.075387
ind_test(as.vector(V99)) #14.92951
ind_test(as.vector(V95))  #31.13947
1-pchisq(ind_test(as.vector(V995)),1)#0.2997319
1-pchisq(ind_test(as.vector(V99)),1)#0.0001116033
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.06094307
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.00
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00

#ES Test
Z2(p=0.005,ES=ES995_hs,L=dax_log_xts[1051:6826],v=V995)#0.568882
Z2(p=0.01,ES=ES99_hs,L=dax_log_xts[1051:6826],v=V99)#0.4667521
Z2(p=0.05,ES=ES95_hs,L=dax_log_xts[1051:6826],v=V95)#0.1867417

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_hs,alpha=0.005)#0.006285753
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_hs,alpha=0.01)#0.00516493
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_hs,alpha=0.05)#0.0005332478

##NV
sp=numeric(0)
VaR95_nv=numeric(0)
VaR99_nv=numeric(0)
VaR995_nv=numeric(0)
ES95_nv=numeric(0)
ES99_nv=numeric(0)
ES995_nv=numeric(0)

for (i in 1:5776){
  sp=dax_log_xts[i:(1049+i)]
  VaR95_nv[i]=mean(sp)+sd(sp)*qnorm(0.95)
  VaR99_nv[i]=mean(sp)+sd(sp)*qnorm(0.99)
  VaR995_nv[i]=mean(sp)+sd(sp)*qnorm(0.995)
  ES95_nv[i]=mean(sp)+sd(sp)*dnorm(qnorm(0.95))/0.05
  ES99_nv[i]=mean(sp)+sd(sp)*dnorm(qnorm(0.99))/0.01
  ES995_nv[i]=mean(sp)+sd(sp)*dnorm(qnorm(0.995))/0.005
}

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_nv,dax_log$date[1051:6826]),col="red")   
lines(xts(VaR99_nv,dax_log$date[1051:6826]),col="blue")    
lines(xts(VaR95_nv,dax_log$date[1051:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_nv<dax_log_xts[1051:6826]) #  88 Ueberschreitungen
sum(VaR99_nv<dax_log_xts[1051:6826]) #  144 Ueberschreitungen
sum(VaR95_nv<dax_log_xts[1051:6826]) #  331 Ueberschreitungen

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,10))  
lines(xts(ES995_nv,dax_log$date[1051:6826]),col="red")   
lines(xts(ES99_nv,dax_log$date[1051:6826]),col="blue")    
lines(xts(ES95_nv,dax_log$date[1051:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_nv<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_nv<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_nv<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #78.46726
uc_test(p=0.01,v=V99)  #91.92006
uc_test(p=0.05,v=V95)  #6.211567
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0
1-pchisq(uc_test(p=0.01,v=V99),1)#0
1-pchisq(uc_test(p=0.05,v=V95),1)#0.013

##Ind.Test 
ind_test(as.vector(V995))#12.57878
ind_test(as.vector(V99)) #25.28661
ind_test(as.vector(V95))  #32.3934
1-pchisq(ind_test(as.vector(V995)),1)#0.00
1-pchisq(ind_test(as.vector(V99)),1)#0.00
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0


####ES Test

Z2(p=0.005,ES=ES995_nv,L=dax_log_xts[1051:6826],v=V995)#2.862474
Z2(p=0.01,ES=ES99_nv,L=dax_log_xts[1051:6826],v=V99)#1.983503
Z2(p=0.05,ES=ES95_nv,L=dax_log_xts[1051:6826],v=V95)#0.3801311

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_nv,alpha=0.005)#0
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_nv,alpha=0.01)#0
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_nv,alpha=0.05)#0


##t
sp=numeric(0)
VaR95_t=numeric(0)
VaR99_t=numeric(0)
VaR995_t=numeric(0)
ES95_t=numeric(0)
ES99_t=numeric(0)
ES995_t=numeric(0)

for (i in 1:5776){
  sp=dax_log_xts[i:(1049+i)]
  tfit=fitdistr(sp,"t")#lower=c(-1, 0.001, 1))
  para=as.numeric(tfit$estimate)
  VaR95_t[i]=para[1]+para[2]*qt(0.95,df=para[3])
  VaR99_t[i]=para[1]+para[2]*qt(0.99,df=para[3])
  VaR995_t[i]=para[1]+para[2]*qt(0.995,df=para[3])
  ES95_t[i]=para[1]+para[2]*(dt(qt(0.95,df=para[3]),df=para[3])/0.05*((para[3]+(qt(0.95,df=para[3]))^2)/(para[3]-1)))
  ES99_t[i]=para[1]+para[2]*(dt(qt(0.99,df=para[3]),df=para[3])/0.01*((para[3]+(qt(0.99,df=para[3]))^2)/(para[3]-1)))
  ES995_t[i]=para[1]+para[2]*(dt(qt(0.995,df=para[3]),df=para[3])/0.005*((para[3]+(qt(0.995,df=para[3]))^2)/(para[3]-1)))
}

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(xts(VaR995_t,dax_log$date[1051:6826]),col="red")   
lines(xts(VaR99_t,dax_log$date[1051:6826]),col="blue")    
lines(xts(VaR95_t,dax_log$date[1051:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_t<dax_log_xts[1051:6826]) #  44 Ueberschreitungen
sum(VaR99_t<dax_log_xts[1051:6826]) #  86 Ueberschreitungen
sum(VaR95_t<dax_log_xts[1051:6826]) #  377 Ueberschreitungen

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
lines(xts(ES995_t,dax_log$date[1051:6826]),col="red")   
lines(xts(ES99_t,dax_log$date[1051:6826]),col="blue")    
lines(xts(ES95_t,dax_log$date[1051:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_t<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_t<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_t<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #6.851362
uc_test(p=0.01,v=V99)  #12.12443
uc_test(p=0.05,v=V95)  #25.97449
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.008857372
1-pchisq(uc_test(p=0.01,v=V99),1)#0
1-pchisq(uc_test(p=0.05,v=V95),1)#0

##Ind.Test 
ind_test(as.vector(V995))#3.943795
ind_test(as.vector(V99)) #16.97729
ind_test(as.vector(V95))  #39.42027
1-pchisq(ind_test(as.vector(V995)),1)#0.04704453
1-pchisq(ind_test(as.vector(V99)),1)#0.00
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.00452753
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.00
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00


###ES TEST

Z2(p=0.005,ES=ES995_t,L=dax_log_xts[1051:6826],v=V995)#0.6126703
Z2(p=0.01,ES=ES99_t,L=dax_log_xts[1051:6826],v=V99)#0.5422024
Z2(p=0.05,ES=ES95_t,L=dax_log_xts[1051:6826],v=V95)#0.3615845

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_t,alpha=0.005)#0.01262904
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_t,alpha=0.01)#0.004443261
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_t,alpha=0.05)#0.00

#skt
sp=numeric(0)
VaR95_skt1=numeric(0)
VaR99_skt1=numeric(0)
VaR995_skt1=numeric(0)
ES95_skt1=numeric(0)
ES99_skt1=numeric(0)
ES995_skt1=numeric(0)

#integrand1=function(x){2*x/(paraskt1[4] + 1/paraskt1[4]) * dt(((x[x >= 0]-paraskt1[1])/paraskt1[2])/paraskt1[4], df=paraskt1[3])/paraskt1[2]}
f1=function(x) qdist("sstd", p = x, mu =paraskt1[1] , sigma = paraskt1[2], skew = 
                       paraskt1[4], shape = paraskt1[3])

mf=function(nu,xi){gamma((nu-1)/2)*sqrt(nu-2)/(sqrt(pi)*gamma(nu/2))*(xi-1/xi)}
sf=function(xi,m){(xi^2+1/xi-1)-m^2}

for (i in 1:5776){
  sp=-dax_log_xts[i:(1049+i)]
  skt1fit=sstdFit(sp)
  paraskt1=as.numeric(skt1fit$estimate)
  VaR95_skt1[i]=qsstd(0.05,mean=paraskt1[1],sd=paraskt1[2],nu=paraskt1[3],xi=paraskt1[4])
  VaR99_skt1[i]=qsstd(0.01,mean=paraskt1[1],sd=paraskt1[2],nu=paraskt1[3],xi=paraskt1[4])
  VaR995_skt1[i]=qsstd(0.005,mean=paraskt1[1],sd=paraskt1[2],nu=paraskt1[3],xi=paraskt1[4])
  ES95_skt1[i]=integrate(f1,lower=0,upper=0.05)$value*20
  ES99_skt1[i]=integrate(f1,lower=0,upper=0.01)$value*100
  ES995_skt1[i]=integrate(f1,lower=0,upper=0.005)$value*200
}
# 1/0.95929623*qt(0.025*(1+0.95929623^2),df=4.00246904)

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(xts(-VaR995_skt1,dax_log$date[1051:6826]),col="red")   
lines(xts(-VaR99_skt1,dax_log$date[1051:6826]),col="blue")    
lines(xts(-VaR95_skt1,dax_log$date[1051:6826]),col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(-VaR995_skt1<dax_log_xts[1051:6826]) #  41 Ueberschreitungen
sum(-VaR99_skt1<dax_log_xts[1051:6826]) #  73 Ueberschreitungen
sum(-VaR95_skt1<dax_log_xts[1051:6826]) #  342 Ueberschreitungen

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
lines(xts(-ES995_skt1,dax_log$date[1051:6826]),col="red")   
lines(xts(-ES99_skt1,dax_log$date[1051:6826]),col="blue")    
lines(xts(-ES95_skt1,dax_log$date[1051:6826]),col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(-VaR995_skt1<dax_log_xts[1051:6826])
sum(V995)
V99=(-VaR99_skt1<dax_log_xts[1051:6826])
sum(V99)
V95=(-VaR95_skt1<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #4.520243
uc_test(p=0.01,v=V99)  #3.748443
uc_test(p=0.05,v=V95)  #9.765674
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.03349607
1-pchisq(uc_test(p=0.01,v=V99),1)#0.05285672
1-pchisq(uc_test(p=0.05,v=V95),1)#0.001778003

##Ind.Test 
ind_test(as.vector(V995))#4.437392
ind_test(as.vector(V99)) #9.216849
ind_test(as.vector(V95))  #30.35321
1-pchisq(ind_test(as.vector(V995)),1)#0.03515992
1-pchisq(ind_test(as.vector(V99)),1)#0.002397979
1-pchisq(ind_test(as.vector(V95)),1)#0.00

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.01134682
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.001529758
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00


####ES TEST

Z2(p=0.005,ES=-ES995_skt1,L=dax_log_xts[1051:6826],v=V995)#0.4232022
Z2(p=0.01,ES=-ES99_skt1,L=dax_log_xts[1051:6826],v=V99)#0.2711923
Z2(p=0.05,ES=-ES95_skt1,L=dax_log_xts[1051:6826],v=V95)#0.2008423

esr_backtest_intercept(-dax_log_xts[1051:6826],e=ES995_skt1,alpha=0.005)#0.09472523
esr_backtest_intercept(-dax_log_xts[1051:6826],e=ES99_skt1,alpha=0.01)#0.08205396
esr_backtest_intercept(-dax_log_xts[1051:6826],e=ES95_skt1,alpha=0.05)#0.006064479

##GARCH mit NV
VaR95_garchnv_z=numeric(0)
VaR99_garchnv_z=numeric(0)
VaR995_garchnv_z=numeric(0)
ES95_garchnv_z=numeric(0)
ES99_garchnv_z=numeric(0)
ES995_garchnv_z=numeric(0)


VaR95_garchnv=numeric(0)
VaR99_garchnv=numeric(0)
VaR995_garchnv=numeric(0)
ES95_garchnv=numeric(0)
ES99_garchnv=numeric(0)
ES995_garchnv=numeric(0)


for (i in (1:5776)){         #es gibt (6826-1050) Vorhersagen. 
  garchnvfit=garchFit(formula=~arma(0,0)+garch(1,1),data=dax_log_xts[i:(1049+i)],cond.dist ="norm")
  zgarchnv=(dax_log_xts[i:(1049+i)]-garchnvfit@fitted)/garchnvfit@sigma.t
  VaR995_garchnv_z=qnorm(0.995,mean=mean(zgarchnv),sd=sd(zgarchnv))
  VaR99_garchnv_z=qnorm(0.99,mean=mean(zgarchnv),sd=sd(zgarchnv))
  VaR95_garchnv_z=qnorm(0.95,mean=mean(zgarchnv),sd=sd(zgarchnv))
  
  VaR995_garchnv[i]=VaR995_garchnv_z*predict(garchnvfit)[1,3]+predict(garchnvfit)[1,1]
  VaR99_garchnv[i]=VaR99_garchnv_z*predict(garchnvfit)[1,3]+predict(garchnvfit)[1,1]
  VaR95_garchnv[i]=VaR95_garchnv_z*predict(garchnvfit)[1,3]+predict(garchnvfit)[1,1]

  ES995_garchnv_z[i]=mean(zgarchnv)+sd(zgarchnv)*dnorm(qnorm(0.995))/0.005
  ES99_garchnv_z[i]=mean(zgarchnv)+sd(zgarchnv)*dnorm(qnorm(0.99))/0.01
  ES95_garchnv_z[i]=mean(zgarchnv)+sd(zgarchnv)*dnorm(qnorm(0.95))/0.05
 
  ES995_garchnv[i]=ES995_garchnv_z[i]*predict(garchnvfit)[1,3]+predict(garchnvfit)[1,1]
  ES99_garchnv[i]=ES99_garchnv_z[i]*predict(garchnvfit)[1,3]+predict(garchnvfit)[1,1]
  ES95_garchnv[i]=ES95_garchnv_z[i]*predict(garchnvfit)[1,3]+predict(garchnvfit)[1,1]
  
}  
VaR995_garchnv_xts=xts(VaR995_garchnv,dax_log$date[1051:6826])
VaR99_garchnv_xts=xts(VaR99_garchnv,dax_log$date[1051:6826])
VaR95_garchnv_xts=xts(VaR95_garchnv,dax_log$date[1051:6826])

ES995_garchnv_xts=xts(ES995_garchnv,dax_log$date[1051:6826])
ES99_garchnv_xts=xts(ES99_garchnv,dax_log$date[1051:6826])
ES95_garchnv_xts=xts(ES95_garchnv,dax_log$date[1051:6826])

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_garchnv_xts,col="red")   
lines(VaR99_garchnv_xts,col="blue")    
lines(VaR95_garchnv_xts,col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_garchnv_xts<dax_log_xts[1051:6826]) #  59 Ueberschreitungen
sum(VaR99_garchnv_xts<dax_log_xts[1051:6826]) #  93 Ueberschreitungen
sum(VaR95_garchnv_xts<dax_log_xts[1051:6826]) #  344 Ueberschreitungen

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
lines(ES995_garchnv_xts,col="red")   
lines(ES99_garchnv_xts,col="blue")    
lines(ES95_garchnv_xts,col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#nur die ersten 1000 Prognosen
plot(dax_log_xts[1051:2050],main="VaR",ylim=c(0,12))  
lines(VaR995_garchnv_xts[1:1000],col="red")   
lines(VaR99_garchnv_xts[1:1000],col="blue")    
lines(VaR95_garchnv_xts[1:1000],col="green")   
legend("topleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

plot(dax_log_xts[1051:2050],main="ES",ylim=c(0,13))  
lines(ES995_garchnv_xts[1:1000],col="red")   
lines(ES99_garchnv_xts[1:1000],col="blue")    
lines(ES95_garchnv_xts[1:1000],col="green") 
legend("topleft",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_garchnv_xts<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_garchnv_xts<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_garchnv_xts<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #24.21593
uc_test(p=0.01,v=V99)  #18.32998
uc_test(p=0.05,v=V95)  #10.49336
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.00
1-pchisq(uc_test(p=0.01,v=V99),1)#0.00
1-pchisq(uc_test(p=0.05,v=V95),1)#0.00

##Ind.Test 
ind_test(as.vector(V995))#0.2234227
ind_test(as.vector(V99)) #0.1578624
ind_test(as.vector(V95))  #0.7092295
1-pchisq(ind_test(as.vector(V995)),1)#0.6364443
1-pchisq(ind_test(as.vector(V99)),1)#0.6911322
1-pchisq(ind_test(as.vector(V95)),1)#0.3996997

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.00
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.00
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.00

##ES TEST

Z2(p=0.005,ES=ES995_garchnv_xts,L=dax_log_xts[1051:6826],v=V995)#1.224931
Z2(p=0.01,ES=ES99_garchnv_xts,L=dax_log_xts[1051:6826],v=V99)#0.7456903
Z2(p=0.05,ES=ES95_garchnv_xts,L=dax_log_xts[1051:6826],v=V95)#0.2678877

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_garchnv_xts,alpha=0.005)#0.0003640444
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_garchnv_xts,alpha=0.01)#0.00
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_garchnv_xts,alpha=0.05)#0.00

##GARCH mit t
VaR95_garcht_z=numeric(0)
VaR99_garcht_z=numeric(0)
VaR995_garcht_z=numeric(0)
ES95_garcht_z=numeric(0)
ES99_garcht_z=numeric(0)
ES995_garcht_z=numeric(0)


VaR95_garcht=numeric(0)
VaR99_garcht=numeric(0)
VaR995_garcht=numeric(0)
ES95_garcht=numeric(0)
ES99_garcht=numeric(0)
ES995_garcht=numeric(0)


for (i in (1:5776)){         #es gibt (6826-1050) Vorhersagen.
  garchtfit=garchFit(formula=~arma(0,0)+garch(1,1),data=dax_log_xts[i:(1049+i)],cond.dist ="std")
  zgarcht=(dax_log_xts[i:(1049+i)]-garchtfit@fitted)/garchtfit@sigma.t
##
  tgarchfit=fitdistr(zgarcht,"t")
  paragarcht=as.numeric(tgarchfit$estimate)
  VaR95_garcht_z[i]=paragarcht[1]+paragarcht[2]*qt(0.95,df=paragarcht[3])
  VaR99_garcht_z[i]=paragarcht[1]+paragarcht[2]*qt(0.99,df=paragarcht[3])
  VaR995_garcht_z[i]=paragarcht[1]+paragarcht[2]*qt(0.995,df=paragarcht[3])
  ES95_garcht_z[i]=paragarcht[1]+paragarcht[2]*(dt(qt(0.95,df=paragarcht[3]),df=paragarcht[3])/0.05*((paragarcht[3]+(qt(0.95,df=paragarcht[3]))^2)/(paragarcht[3]-1)))
  ES99_garcht_z[i]=paragarcht[1]+paragarcht[2]*(dt(qt(0.99,df=paragarcht[3]),df=paragarcht[3])/0.01*((paragarcht[3]+(qt(0.99,df=paragarcht[3]))^2)/(paragarcht[3]-1)))
  ES995_garcht_z[i]=paragarcht[1]+paragarcht[2]*(dt(qt(0.995,df=paragarcht[3]),df=paragarcht[3])/0.005*((paragarcht[3]+(qt(0.995,df=paragarcht[3]))^2)/(paragarcht[3]-1)))  
  
  VaR995_garcht[i]=VaR995_garcht_z[i]*predict(garchtfit)[1,3]+predict(garchtfit)[1,1]
  VaR99_garcht[i]=VaR99_garcht_z[i]*predict(garchtfit)[1,3]+predict(garchtfit)[1,1]
  VaR95_garcht[i]=VaR95_garcht_z[i]*predict(garchtfit)[1,3]+predict(garchtfit)[1,1]
  
  
  ES995_garcht[i]=ES995_garcht_z[i]*predict(garchtfit)[1,3]+predict(garchtfit)[1,1]
  ES99_garcht[i]=ES99_garcht_z[i]*predict(garchtfit)[1,3]+predict(garchtfit)[1,1]
  ES95_garcht[i]=ES95_garcht_z[i]*predict(garchtfit)[1,3]+predict(garchtfit)[1,1]
  
}  
VaR995_garcht_xts=xts(VaR995_garcht,dax_log$date[1051:6826])
VaR99_garcht_xts=xts(VaR99_garcht,dax_log$date[1051:6826])
VaR95_garcht_xts=xts(VaR95_garcht,dax_log$date[1051:6826])

ES995_garcht_xts=xts(ES995_garcht,dax_log$date[1051:6826])
ES99_garcht_xts=xts(ES99_garcht,dax_log$date[1051:6826])
ES95_garcht_xts=xts(ES95_garcht,dax_log$date[1051:6826])

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
lines(VaR995_garcht_xts,col="red")   
lines(VaR99_garcht_xts,col="blue")    
lines(VaR95_garcht_xts,col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_garcht_xts<dax_log_xts[1051:6826]) #  36 Ueberschreitungen
sum(VaR99_garcht_xts<dax_log_xts[1051:6826]) #  72 Ueberschreitungen
sum(VaR95_garcht_xts<dax_log_xts[1051:6826]) #  366 Ueberschreitungen

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
lines(ES995_garcht_xts,col="red")   
lines(ES99_garcht_xts,col="blue")    
lines(ES95_garcht_xts,col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#nur die ersten 1000 Prognosen
plot(dax_log_xts[1051:2050],main="VaR",ylim=c(0,12))  
lines(VaR995_garcht_xts[1:1000],col="red")   
lines(VaR99_garcht_xts[1:1000],col="blue")    
lines(VaR95_garcht_xts[1:1000],col="green")   
legend("topleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

plot(dax_log_xts[1051:2050],main="ES",ylim=c(0,13))  
lines(ES995_garcht_xts[1:1000],col="red")   
lines(ES99_garcht_xts[1:1000],col="blue")    
lines(ES95_garcht_xts[1:1000],col="green") 
legend("topleft",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_garcht_xts<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_garcht_xts<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_garcht_xts<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #1.635437
uc_test(p=0.01,v=V99)  #3.288717
uc_test(p=0.05,v=V95)  #20.10128
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.2009526
1-pchisq(uc_test(p=0.01,v=V99),1)#0.06975751
1-pchisq(uc_test(p=0.05,v=V95),1)#0.00

##Ind.Test 
ind_test(as.vector(V995))#1.471322
ind_test(as.vector(V99)) #0.01154386
ind_test(as.vector(V95))  #0.5234736
1-pchisq(ind_test(as.vector(V995)),1)#0.2251373
1-pchisq(ind_test(as.vector(V99)),1)#0.9144381
1-pchisq(ind_test(as.vector(V95)),1)#0.4693637

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.2115319
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.1920249
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.000

##ES TEST

Z2(p=0.005,ES=ES995_garcht_xts,L=dax_log_xts[1051:6826],v=V995)#0.2841605
Z2(p=0.01,ES=ES99_garcht_xts,L=dax_log_xts[1051:6826],v=V99)#0.2586691
Z2(p=0.05,ES=ES95_garcht_xts,L=dax_log_xts[1051:6826],v=V95)#0.2768927

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_garcht_xts,alpha=0.005)#0.07982554
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_garcht_xts,alpha=0.01)#0.09269161
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_garcht_xts,alpha=0.05)#0.0006222125

###Garch mit skt ##dauert 2 Stunden

#VaR95_garchskt_z=numeric(0)
#VaR99_garchskt_z=numeric(0)
#VaR995_garchskt_z=numeric(0)
#ES95_garchskt_z=numeric(0)
#ES99_garchskt_z=numeric(0)
#ES995_garchskt_z=numeric(0)


#VaR95_garchskt=numeric(0)
#VaR99_garchskt=numeric(0)
#VaR995_garchskt=numeric(0)
#ES95_garchskt=numeric(0)
#ES99_garchskt=numeric(0)
#ES995_garchskt=numeric(0)

#integrandgarchskt=function(x){2*x/(paragarchskt[4] + 1/paragarchskt[4]) * dt(x[x >= 0]/paragarchskt[4], df=paragarchskt[3])}


#for (i in (1:5776)){         #es gibt (6826-1050) Vorhersagen.  
#  garchsktfit=garchFit(formula=~arma(0,0)+garch(1,1),data=dax_log_xts[i:(1049+i)],cond.dist ="sstd")
#  zgarchskt=(dax_log_xts[i:(1049+i)]-garchsktfit@fitted)/garchsktfit@sigma.t
  ##
#  sktgarchfit=sstdFit(zgarchskt)
#  paragarchskt=as.numeric(sktgarchfit$estimate)
#  VaR95_garchskt_z[i]=qskt(0.95,df=paragarchskt[3],gamma=paragarchskt[4])
#  VaR99_garchskt_z[i]=qskt(0.99,df=paragarchskt[3],gamma=paragarchskt[4])
#  VaR995_garchskt_z[i]=qskt(0.995,df=paragarchskt[3],gamma=paragarchskt[4])
#  ES95_garchskt_z[i]=integrate(integrandgarchskt,lower=qskt(0.95,df=paragarchskt[3],gamma=paragarchskt[4]),upper=Inf)$value*20
#  ES99_garchskt_z[i]=integrate(integrandgarchskt,lower=qskt(0.99,df=paragarchskt[3],gamma=paragarchskt[4]),upper=Inf)$value*100
#  ES995_garchskt_z[i]=integrate(integrandgarchskt,lower=qskt(0.995,df=paragarchskt[3],gamma=paragarchskt[4]),upper=Inf)$value*200
  
#  VaR995_garchskt[i]=VaR995_garchskt_z[i]*predict(garchsktfit)[1,3]+predict(garchsktfit)[1,1]
#  VaR99_garchskt[i]=VaR99_garchskt_z[i]*predict(garchsktfit)[1,3]+predict(garchsktfit)[1,1]
#  VaR95_garchskt[i]=VaR95_garchskt_z[i]*predict(garchsktfit)[1,3]+predict(garchsktfit)[1,1]
  
#  ES995_garchskt[i]=ES995_garchskt_z[i]*predict(garchsktfit)[1,3]+predict(garchsktfit)[1,1]
#  ES99_garchskt[i]=ES99_garchskt_z[i]*predict(garchsktfit)[1,3]+predict(garchsktfit)[1,1]
#  ES95_garchskt[i]=ES95_garchskt_z[i]*predict(garchsktfit)[1,3]+predict(garchsktfit)[1,1]
  
#}  
#VaR995_garchskt_xts=xts(VaR995_garchskt,dax_log$date[1051:6826])
#VaR99_garchskt_xts=xts(VaR99_garchskt,dax_log$date[1051:6826])
#VaR95_garchskt_xts=xts(VaR95_garchskt,dax_log$date[1051:6826])

#ES995_garchskt_xts=xts(ES995_garchskt,dax_log$date[1051:6826])
#ES99_garchskt_xts=xts(ES99_garchskt,dax_log$date[1051:6826])
#ES95_garchskt_xts=xts(ES95_garchskt,dax_log$date[1051:6826])

#plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,10))  
#lines(VaR995_garchskt_xts,col="red")   
#lines(VaR99_garchskt_xts,col="blue")    
#lines(VaR95_garchskt_xts,col="green")   
#legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
#sum(VaR995_garchskt_xts<dax_log_xts[1051:6826]) #  10 Ueberschreitungen
#sum(VaR99_garchskt_xts<dax_log_xts[1051:6826]) #  15 Ueberschreitungen
#sum(VaR95_garchskt_xts<dax_log_xts[1051:6826]) #  132 Ueberschreitungen

#plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,12))  
#lines(ES995_garchskt_xts,col="red")   
#lines(ES99_garchskt_xts,col="blue")    
#lines(ES95_garchskt_xts,col="green") 
#legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))#

#nur die ersten 1000 Prognosen
#plot(dax_log_xts[1051:2050],main="VaR",ylim=c(0,12))  
#lines(VaR995_garchskt_xts[1:1000],col="red")   
#lines(VaR99_garchskt_xts[1:1000],col="blue")    
#lines(VaR95_garchskt_xts[1:1000],col="green")   
#legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#plot(dax_log_xts[1051:2050],main="ES",ylim=c(0,13))  
#lines(ES995_garchskt_xts[1:1000],col="red")   
#lines(ES99_garchskt_xts[1:1000],col="blue")    
#lines(ES95_garchskt_xts[1:1000],col="green") 
#legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
#V995=(VaR995_garchskt_xts<dax_log_xts[1051:6826])
#sum(V995)
#V99=(VaR99_garchskt_xts<dax_log_xts[1051:6826])
#sum(V99)
#V95=(VaR95_garchskt_xts<dax_log_xts[1051:6826])
#sum(V95)

#uc_test(p=0.005,v=V995) #8.406399
#uc_test(p=0.01,v=V99)  #26.25402
#uc_test(p=0.05,v=V95)  #44.06991
#p-wert
#1-pchisq(uc_test(p=0.005,v=V995),1)#0.003739026
#1-pchisq(uc_test(p=0.01,v=V99),1)#0.00
#1-pchisq(uc_test(p=0.05,v=V95),1)#0.00

##Ind.Test 
#ind_test(as.vector(V995))#0.04530015
#ind_test(as.vector(V99)) #4.182913
#ind_test(as.vector(V95))  #0.2573659
#1-pchisq(ind_test(as.vector(V995)),1)#0.8314531
#1-pchisq(ind_test(as.vector(V99)),1)#0.04083347
#1-pchisq(ind_test(as.vector(V95)),1)#0.6119356

##CC Test
#1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.01461291
#1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.00
#1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.000


##ES TEST
#Z2(p=0.005,ES=ES995_garchskt_xts,L=dax_log_xts[1051:6826],v=V995)#0.4106186
#Z2(p=0.01,ES=ES99_garchskt_xts,L=dax_log_xts[1051:6826],v=V99)#0.1826519
#Z2(p=0.05,ES=ES95_garchskt_xts,L=dax_log_xts[1051:6826],v=V95)#0.2890399

#esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_garchskt_xts,alpha=0.005)#0.06637217
#esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_garchskt_xts,alpha=0.01)#0.08002596
#esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_garchskt_xts,alpha=0.05)#0.002496068

###GARCH-Skewed t (neu)
model=ugarchspec( variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "sstd" ) 
modelroll=ugarchroll(model,data=-dax_log_xts,n.ahead = 1,forecast.length = 5776,window.size = 1050,refit.every = 1,calculate.VaR = TRUE,VaR.alpha = c(0.005,0.01,0.05),keep.coef=TRUE,cluster=NULL)
#modelfit
modelroll@forecast
Hit=modelroll@forecast$VaR[,"realized"]<modelroll@forecast$VaR[,"alpha(0%)"]
sum(Hit)
modelroll@forecast

VaR995_garchsktnew_xts=xts(-modelroll@forecast$VaR[,"alpha(0%)"],dax_log$date[1051:6826])
VaR99_garchsktnew_xts=xts(-modelroll@forecast$VaR[,"alpha(1%)"],dax_log$date[1051:6826])
VaR95_garchsktnew_xts=xts(-modelroll@forecast$VaR[,"alpha(5%)"],dax_log$date[1051:6826])

plot(dax_log_xts[1051:6826],main="VaR",ylim=c(0,12))  
lines(VaR995_garchsktnew_xts,col="red")   
lines(VaR99_garchsktnew_xts,col="blue")    
lines(VaR95_garchsktnew_xts,col="green")   
legend("topright",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_garchsktnew_xts<dax_log_xts[1051:6826]) #  26 Ueberschreitungen
sum(VaR99_garchsktnew_xts<dax_log_xts[1051:6826]) #  61 Ueberschreitungen
sum(VaR95_garchsktnew_xts<dax_log_xts[1051:6826]) #  353 Ueberschreitungen

#nur die ersten 1000 Prognosen
plot(dax_log_xts[1051:2050],main="VaR",ylim=c(0,12))  
lines(VaR995_garchsktnew_xts[1:1000],col="red")   
lines(VaR99_garchsktnew_xts[1:1000],col="blue")    
lines(VaR95_garchsktnew_xts[1:1000],col="green")   
legend("topleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

##ES
ES95_garchsktnew=numeric(0)
ES99_garchsktnew=numeric(0)
ES995_garchsktnew=numeric(0)

df1_var <- as.data.frame(modelroll, which = "density") 
f = function(x, skew, shape) qdist("sstd", p = x, mu = 0, sigma = 1, skew = 
                                     skew, shape = shape) 
ES95_garchsktnew = -apply(df1_var, 1, function(x) x['Mu'] + x['Sigma'] * integrate(f,0, 0.05, skew = x['Skew'], shape = x['Shape'])$value*20) 
ES99_garchsktnew = -apply(df1_var, 1, function(x) x['Mu'] + x['Sigma'] * integrate(f,0, 0.01, skew = x['Skew'], shape = x['Shape'])$value*100) 
ES995_garchsktnew = -apply(df1_var, 1, function(x) x['Mu'] + x['Sigma'] * integrate(f,0, 0.005, skew = x['Skew'], shape = x['Shape'])$value*200) 

#for (i in (1:5776)){
#ES95_garchsktnew[i]=integrate(integrandgarchskt,lower=qskt(0.95,df=modelroll@forecast$density$Shape[i],gamma=modelroll@forecast$density$Skew[i]),upper=Inf)$value*20
#ES99_garchsktnew[i]=integrate(integrandgarchskt,lower=qskt(0.99,df=modelroll@forecast$density$Shape[i],gamma=modelroll@forecast$density$Skew[i]),upper=Inf)$value*100
#ES995_garchsktnew[i]=integrate(integrandgarchskt,lower=qskt(0.995,df=modelroll@forecast$density$Shape[i],gamma=modelroll@forecast$density$Skew[i]),upper=Inf)$value*200
#}
ES95_garchsktnew_xts=xts(ES95_garchsktnew,dax_log$date[1051:6826])
ES99_garchsktnew_xts=xts(ES99_garchsktnew,dax_log$date[1051:6826])
ES995_garchsktnew_xts=xts(ES995_garchsktnew,dax_log$date[1051:6826])

plot(dax_log_xts[1051:6826],main="ES",ylim=c(0,15))  
lines(ES995_garchsktnew_xts,col="red")   
lines(ES99_garchsktnew_xts,col="blue")    
lines(ES95_garchsktnew_xts,col="green") 
legend("topright",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#nur die ersten 1000
plot(dax_log_xts[1051:2050],main="ES",ylim=c(0,15))  
lines(ES995_garchsktnew_xts[1:1000],col="red")   
lines(ES99_garchsktnew_xts[1:1000],col="blue")    
lines(ES95_garchsktnew_xts[1:1000],col="green") 
legend("topleft",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

#U.C Test
V995=(VaR995_garchsktnew_xts<dax_log_xts[1051:6826])
sum(V995)
V99=(VaR99_garchsktnew_xts<dax_log_xts[1051:6826])
sum(V99)
V95=(VaR95_garchsktnew_xts<dax_log_xts[1051:6826])
sum(V95)

uc_test(p=0.005,v=V995) #0.2986986
uc_test(p=0.01,v=V99)  #0.1802752
uc_test(p=0.05,v=V95)  #14.07204
#p-wert
1-pchisq(uc_test(p=0.005,v=V995),1)#0.5846994
1-pchisq(uc_test(p=0.01,v=V99),1)#0.6711368
1-pchisq(uc_test(p=0.05,v=V95),1)#0.0001759398

##Ind.Test 
ind_test(as.vector(V995))#2.585386
ind_test(as.vector(V99)) #0.1719721
ind_test(as.vector(V95))  #0.1336985
1-pchisq(ind_test(as.vector(V995)),1)#0.1078541
1-pchisq(ind_test(as.vector(V99)),1)#0.678365
1-pchisq(ind_test(as.vector(V95)),1)#0.7146277

##CC Test
1-pchisq(uc_test(p=0.005,v=V995)+ind_test(as.vector(V995)),2)#0.2364443
1-pchisq(uc_test(p=0.01,v=V99)+ind_test(as.vector(V99)),2)#0.8385143
1-pchisq(uc_test(p=0.05,v=V95)+ind_test(as.vector(V95)),2)#0.0008227406


###ES TEST

Z2(p=0.005,ES=ES995_garchsktnew_xts,L=dax_log_xts[1051:6826],v=V995)#-0.09013499
Z2(p=0.01,ES=ES99_garchsktnew_xts,L=dax_log_xts[1051:6826],v=V99)#0.03357583
Z2(p=0.05,ES=ES95_garchsktnew_xts,L=dax_log_xts[1051:6826],v=V95)#0.1965791

esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES995_garchsktnew_xts,alpha=0.005)#0.4586323
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES99_garchsktnew_xts,alpha=0.01)#0.5619296
esr_backtest_intercept(-dax_log_xts[1051:6826],e=-ES95_garchsktnew_xts,alpha=0.05)#0.07076329

###Multidays VaR
length(dax_xts)
length(dax_log_xts)
-100*log(dax$Close[11]/dax$Close[1])
sum(dax_log_xts[1:10])
###Cheng et al. 2011 S.1160
##h=10,T=6000,s=1:599,k=1:10 (in Literatur ist k=1:9)
dax_log_xts6000=dax_log_xts[1:6000]
#Subsample 1 k=1
s1=numeric(0)
sum(dax_log_xts6000[1:10])
for (i in (1:599)){
  s1[i]=sum(dax_log_xts6000[(10*i-9):(10*i)])
}
length(s1)
head(s1)
#Subsample 2 k=2
s2=numeric(0)
for (i in (1:599)){
  s2[i]=sum(dax_log_xts6000[(10*i-8):(10*i+1)])
}
length(s2)
head(s2)
#Subsample 3 k=3
s3=numeric(0)
for (i in (1:599)){
  s3[i]=sum(dax_log_xts6000[(10*i-7):(10*i+2)])
}
length(s3)
head(s3)
#Subsample 4 k=4
s4=numeric(0)
for (i in (1:599)){
  s4[i]=sum(dax_log_xts6000[(10*i-6):(10*i+3)])
}
length(s4)
head(s4)
#Subsample 5 k=5
s5=numeric(0)
for (i in (1:599)){
  s5[i]=sum(dax_log_xts6000[(10*i-5):(10*i+4)])
}
length(s5)
head(s5)
#Subsample 6 k=6
s6=numeric(0)
for (i in (1:599)){
  s6[i]=sum(dax_log_xts6000[(10*i-4):(10*i+5)])
}
length(s6)
head(s6)
#Subsample 7 k=7
s7=numeric(0)
for (i in (1:599)){
  s7[i]=sum(dax_log_xts6000[(10*i-3):(10*i+6)])
}
length(s7)
head(s7)
#Subsample 8 k=8
s8=numeric(0)
for (i in (1:599)){
  s8[i]=sum(dax_log_xts6000[(10*i-2):(10*i+7)])
}
length(s8)
head(s8)
#Subsample 9 k=9
s9=numeric(0)
for (i in (1:599)){
  s9[i]=sum(dax_log_xts6000[(10*i-1):(10*i+8)])
}
length(s9)
head(s9)
#Subsample 10 k=10
s10=numeric(0)
for (i in (1:599)){
  s10[i]=sum(dax_log_xts6000[(10*i):(10*i+9)])
}
length(s10)
head(s10)

#VaR(1) taeglich h=1
quantile(dax_log_xts6000,0.95) #2.29 10^(0.5)*VaR(1)=7.23
quantile(dax_log_xts6000,0.99) #4.31 10^(0.5)*VaR(1)=13.64
quantile(dax_log_xts6000,0.995)#5.15 10^(0.5)*VaR(1)=16.28

###VaR(h)=1/h*Sum_(k=1)^(k=10)VaR_k(h)
##alpha=0.05
(quantile(s1,0.95)+quantile(s2,0.95)+quantile(s3,0.95)+quantile(s4,0.95)+quantile(s5,0.95)+
  quantile(s6,0.95)+quantile(s7,0.95)+quantile(s8,0.95)+quantile(s9,0.95)+quantile(s10,0.95))/10
##6.75  S.1160 zeta(h)=7.11
##alpha=0.01
(quantile(s1,0.99)+quantile(s2,0.99)+quantile(s3,0.99)+quantile(s4,0.99)+quantile(s5,0.99)+
  quantile(s6,0.99)+quantile(s7,0.99)+quantile(s8,0.99)+quantile(s9,0.99)+quantile(s10,0.99))/10
##13.34  zeta(h)=2.25
##alpha=0.005
(quantile(s1,0.995)+quantile(s2,0.995)+quantile(s3,0.995)+quantile(s4,0.995)+quantile(s5,0.995)+
  quantile(s6,0.995)+quantile(s7,0.995)+quantile(s8,0.995)+quantile(s9,0.995)+quantile(s10,0.995))/10
##16.16 zeta(h)=0.74

##h=30,T=6000,s=1:199,k=1:30 (in Literatur ist k=1:29)
#Subsample 1 k=1
s30=matrix(0,nrow=199,ncol=30)
#for (i in (1:199)){
#  s1_30[i]=sum(dax_log_xts6000[(30*i-29):(30*i)])
#}
#length(s1_30)
#head(s1_30)

for (k in (1:30)){
  for (i in (1:199)){
    s30[i,k]=sum(dax_log_xts6000[(30*i-30+k):(30*i+k-1)])
  }
}
s30

s30[,1]

#VaR(1) taeglich h=1
quantile(dax_log_xts6000,0.95) #2.29 30^(0.5)*VaR(1)=12.53
quantile(dax_log_xts6000,0.99) #4.31 30^(0.5)*VaR(1)=23.62
quantile(dax_log_xts6000,0.995)#5.15 30^(0.5)*VaR(1)=28.19

###VaR(h)=1/h*Sum_(k=1)^(k=10)VaR_k(h)
##alpha=0.05
(quantile(s30[,1],0.95)+quantile(s30[,2],0.95)+quantile(s30[,3],0.95)+quantile(s30[,4],0.95)+
quantile(s30[,5],0.95)+quantile(s30[,6],0.95)+quantile(s30[,7],0.95)+quantile(s30[,8],0.95)+
quantile(s30[,9],0.95)+quantile(s30[,10],0.95)+quantile(s30[,11],0.95)+quantile(s30[,12],0.95)+
quantile(s30[,13],0.95)+quantile(s30[,14],0.95)+quantile(s30[,15],0.95)+quantile(s30[,16],0.95)+
quantile(s30[,17],0.95)+quantile(s30[,18],0.95)+quantile(s30[,19],0.95)+quantile(s30[,20],0.95)+
quantile(s30[,21],0.95)+quantile(s30[,22],0.95)+quantile(s30[,23],0.95)+quantile(s30[,24],0.95)+
quantile(s30[,25],0.95)+quantile(s30[,26],0.95)+quantile(s30[,27],0.95)+quantile(s30[,28],0.95)+
quantile(s30[,29],0.95)+quantile(s30[,30],0.95))/30
##12.13 S.1160 zeta(h)=100*((12.53/12.13)-1)=3.30
##alpha=0.01
(quantile(s30[,1],0.99)+quantile(s30[,2],0.99)+quantile(s30[,3],0.99)+quantile(s30[,4],0.99)+
  quantile(s30[,5],0.99)+quantile(s30[,6],0.99)+quantile(s30[,7],0.99)+quantile(s30[,8],0.99)+
  quantile(s30[,9],0.99)+quantile(s30[,10],0.99)+quantile(s30[,11],0.99)+quantile(s30[,12],0.99)+
  quantile(s30[,13],0.99)+quantile(s30[,14],0.99)+quantile(s30[,15],0.99)+quantile(s30[,16],0.99)+
  quantile(s30[,17],0.99)+quantile(s30[,18],0.99)+quantile(s30[,19],0.99)+quantile(s30[,20],0.99)+
  quantile(s30[,21],0.99)+quantile(s30[,22],0.99)+quantile(s30[,23],0.99)+quantile(s30[,24],0.99)+
  quantile(s30[,25],0.99)+quantile(s30[,26],0.99)+quantile(s30[,27],0.99)+quantile(s30[,28],0.99)+
  quantile(s30[,29],0.99)+quantile(s30[,30],0.99))/30
#22.87 zeta(h)=3.28
##alpha=0.005
(quantile(s30[,1],0.995)+quantile(s30[,2],0.995)+quantile(s30[,3],0.995)+quantile(s30[,4],0.995)+
  quantile(s30[,5],0.995)+quantile(s30[,6],0.995)+quantile(s30[,7],0.995)+quantile(s30[,8],0.995)+
  quantile(s30[,9],0.995)+quantile(s30[,10],0.995)+quantile(s30[,11],0.995)+quantile(s30[,12],0.995)+
  quantile(s30[,13],0.995)+quantile(s30[,14],0.995)+quantile(s30[,15],0.995)+quantile(s30[,16],0.995)+
  quantile(s30[,17],0.995)+quantile(s30[,18],0.995)+quantile(s30[,19],0.995)+quantile(s30[,20],0.995)+
  quantile(s30[,21],0.995)+quantile(s30[,22],0.995)+quantile(s30[,23],0.995)+quantile(s30[,24],0.995)+
  quantile(s30[,25],0.995)+quantile(s30[,26],0.995)+quantile(s30[,27],0.995)+quantile(s30[,28],0.995)+
  quantile(s30[,29],0.995)+quantile(s30[,30],0.995))/30
#25.55 zeta(h)=10.33

##Variance-Ratio-Test (Cheng et al. 2011 S.1164)

##k=1:10,h=10,
VR10=((1-1/10)*acf(dax_log_xts6000)$acf[2]+(1-2/10)*acf(dax_log_xts6000)$acf[3]+
  (1-3/10)*acf(dax_log_xts6000)$acf[4]+(1-4/10)*acf(dax_log_xts6000)$acf[5]+
  (1-5/10)*acf(dax_log_xts6000)$acf[6]+(1-6/10)*acf(dax_log_xts6000)$acf[7]+
  (1-7/10)*acf(dax_log_xts6000)$acf[8]+(1-8/10)*acf(dax_log_xts6000)$acf[9]+
  (1-9/10)*acf(dax_log_xts6000)$acf[10])*2+1
VR10 #0.8929

##k=1:30,h=30
VR30=((1-1/30)*acf(dax_log_xts6000)$acf[2]+(1-2/30)*acf(dax_log_xts6000)$acf[3]+
        (1-3/30)*acf(dax_log_xts6000)$acf[4]+(1-4/30)*acf(dax_log_xts6000)$acf[5]+
        (1-5/30)*acf(dax_log_xts6000)$acf[6]+(1-6/30)*acf(dax_log_xts6000)$acf[7]+
        (1-7/30)*acf(dax_log_xts6000)$acf[8]+(1-8/30)*acf(dax_log_xts6000)$acf[9]+
        (1-9/30)*acf(dax_log_xts6000)$acf[10]+(1-10/30)*acf(dax_log_xts6000)$acf[11]+
  (1-11/30)*acf(dax_log_xts6000)$acf[12]+(1-12/30)*acf(dax_log_xts6000)$acf[13]+
  (1-13/30)*acf(dax_log_xts6000)$acf[14]+(1-14/30)*acf(dax_log_xts6000)$acf[15]+
  (1-15/30)*acf(dax_log_xts6000)$acf[16]+(1-16/30)*acf(dax_log_xts6000)$acf[17]+
  (1-17/30)*acf(dax_log_xts6000)$acf[18]+(1-18/30)*acf(dax_log_xts6000)$acf[19]+
  (1-19/30)*acf(dax_log_xts6000)$acf[20]+(1-20/30)*acf(dax_log_xts6000)$acf[21]+
  (1-21/30)*acf(dax_log_xts6000)$acf[22]+(1-22/30)*acf(dax_log_xts6000)$acf[23]+
  (1-23/30)*acf(dax_log_xts6000)$acf[24]+(1-24/30)*acf(dax_log_xts6000)$acf[25]+
  (1-25/30)*acf(dax_log_xts6000)$acf[26]+(1-26/30)*acf(dax_log_xts6000)$acf[27]+
  (1-27/30)*acf(dax_log_xts6000)$acf[28]+(1-28/30)*acf(dax_log_xts6000)$acf[29]+
  (1-29/30)*acf(dax_log_xts6000)$acf[30])*2+1
VR30 #0.8968

###delta_k (S.1164) k=1:30
#mudach=mean(dax_log_xts6000)
#delta_k=numeric(0)
#kj=numeric(0)
#for (k in (1:30)){
#  for (i in ((k+1):6000)){
#   kj[i-k]=(-dax_log$logreturn[i]-mudach)^2*(-dax_log$logreturn[i-k]-mudach)^2
#   delta_k[k]=sum(kj)*6000*(var(dax_log_xts6000)*5999)^(-2)
#    }
#}
#delta_k

#vv=numeric(0)
#for (s in (1:6000)){
#vv[s]=(-dax_log$logreturn[s]-mudach)^2  
#}
#sum(vv)

#vartheta10=((1-1/10)^2*delta_k[1]+(1-2/10)^2*delta_k[2]+
#              (1-3/10)^2*delta_k[3]+(1-4/10)^2*delta_k[4]+
#              (1-5/10)^2*delta_k[5]+(1-6/10)^2*delta_k[6]+
#              (1-7/10)^2*delta_k[7]+(1-8/10)^2*delta_k[8]+
#              (1-9/10)^2*delta_k[9])*4
#vartheta10

#vartheta30=((1-1/30)^2*delta_k[1]+(1-2/30)^2*delta_k[2]+
#            (1-3/30)^2*delta_k[3]+(1-4/30)^2*delta_k[4]+
#            (1-5/30)^2*delta_k[5]+(1-6/30)^2*delta_k[6]+
#            (1-7/30)^2*delta_k[7]+(1-8/30)^2*delta_k[8]+
#            (1-9/30)^2*delta_k[9]+(1-10/30)^2*delta_k[10]+
#            (1-11/30)^2*delta_k[11]+(1-12/30)^2*delta_k[12]+
#            (1-13/30)^2*delta_k[13]+(1-14/30)^2*delta_k[14]+
#            (1-15/30)^2*delta_k[15]+(1-16/30)^2*delta_k[16]+
#            (1-17/30)^2*delta_k[17]+(1-18/30)^2*delta_k[18]+
#            (1-19/30)^2*delta_k[19]+(1-20/30)^2*delta_k[20]+
#            (1-21/30)^2*delta_k[21]+(1-22/30)^2*delta_k[22]+
#            (1-23/30)^2*delta_k[23]+(1-24/30)^2*delta_k[24]+
#            (1-25/30)^2*delta_k[25]+(1-26/30)^2*delta_k[26]+
#            (1-27/30)^2*delta_k[27]+(1-28/30)^2*delta_k[28]+
#            (1-29/30)^2*delta_k[29])*4
#vartheta30

#psi10=sqrt(10)*(VR10-1)/sqrt(vartheta10)
#psi10  #-0.06

#psi30=sqrt(30)*(VR30-1)/sqrt(vartheta30)
#psi30  #-0.0

#Lo-Mackinlay Test (Cheng et al. 2011 S.1164 und Lo, Mackinlay 1988 S.50) 
Lo.Mac(dax_log_xts6000,kvec=c(10,30))  
2*pnorm(-1.5576643)   #h=10,p=0.12
2*pnorm(-0.8459866)   #h=30,p=0.40

#subsample-based Test (Cheng et al. 2011 S.1164)
dax_log_xts6000
daxmat10=matrix(0,nrow=10,ncol=5991)
for(i in (1:5991)){
  daxmat10[,i]=dax_log_xts6000[i:(i+9)]
}
#dax_log_xts6000[1:20]
#daxmat10[,1:2]
ccc=matrix(0,nrow=600,ncol=5000)
resam=matrix(0,nrow=6000,ncol=5000)
var500095=numeric(0)
var500099=numeric(0)
var5000995=numeric(0)
varh500095=numeric(0)
varh500099=numeric(0)
varh5000995=numeric(0)
for (i in 1:5000){
  ccc[,i]=sample(c(1:5991),600,replace = FALSE)
  resam[,i]=daxmat10[,ccc[,i]]
  var500095[i]=quantile(resam[,i],0.95)
  var500099[i]=quantile(resam[,i],0.99)
  var5000995[i]=quantile(resam[,i],0.995)
}

#h=30
sh30=array(0,c(199,30,5000))  #199=nrow, 30=ncol, 5000=n Wiederholung
for (q in 1:5000){
  for (k in (1:30)){
    for (i in (1:199)){
      sh30[i,k,q]=sum(resam[,q][(30*i-30+k):(30*i+k-1)])
    }
  }
}

for (q in 1:5000){
 varh500095[q]= (quantile(sh30[,1,q],0.95)+quantile(sh30[,2,q],0.95)+quantile(sh30[,3,q],0.95)+quantile(sh30[,4,q],0.95)+
     quantile(sh30[,5,q],0.95)+quantile(sh30[,6,q],0.95)+quantile(sh30[,7,q],0.95)+quantile(sh30[,8,q],0.95)+
     quantile(sh30[,9,q],0.95)+quantile(sh30[,10,q],0.95)+quantile(sh30[,11,q],0.95)+quantile(sh30[,12,q],0.95)+
     quantile(sh30[,13,q],0.95)+quantile(sh30[,14,q],0.95)+quantile(sh30[,15,q],0.95)+quantile(sh30[,16,q],0.95)+
     quantile(sh30[,17,q],0.95)+quantile(sh30[,18,q],0.95)+quantile(sh30[,19,q],0.95)+quantile(sh30[,20,q],0.95)+
     quantile(sh30[,21,q],0.95)+quantile(sh30[,22,q],0.95)+quantile(sh30[,23,q],0.95)+quantile(sh30[,24,q],0.95)+
     quantile(sh30[,25,q],0.95)+quantile(sh30[,26,q],0.95)+quantile(sh30[,27,q],0.95)+quantile(sh30[,28,q],0.95)+
     quantile(sh30[,29,q],0.95)+quantile(sh30[,30,q],0.95))/30
 varh500099[q]= (quantile(sh30[,1,q],0.99)+quantile(sh30[,2,q],0.99)+quantile(sh30[,3,q],0.99)+quantile(sh30[,4,q],0.99)+
                   quantile(sh30[,5,q],0.99)+quantile(sh30[,6,q],0.99)+quantile(sh30[,7,q],0.99)+quantile(sh30[,8,q],0.99)+
                   quantile(sh30[,9,q],0.99)+quantile(sh30[,10,q],0.99)+quantile(sh30[,11,q],0.99)+quantile(sh30[,12,q],0.99)+
                   quantile(sh30[,13,q],0.99)+quantile(sh30[,14,q],0.99)+quantile(sh30[,15,q],0.99)+quantile(sh30[,16,q],0.99)+
                   quantile(sh30[,17,q],0.99)+quantile(sh30[,18,q],0.99)+quantile(sh30[,19,q],0.99)+quantile(sh30[,20,q],0.99)+
                   quantile(sh30[,21,q],0.99)+quantile(sh30[,22,q],0.99)+quantile(sh30[,23,q],0.99)+quantile(sh30[,24,q],0.99)+
                   quantile(sh30[,25,q],0.99)+quantile(sh30[,26,q],0.99)+quantile(sh30[,27,q],0.99)+quantile(sh30[,28,q],0.99)+
                   quantile(sh30[,29,q],0.99)+quantile(sh30[,30,q],0.99))/30
 varh5000995[q]= (quantile(sh30[,1,q],0.995)+quantile(sh30[,2,q],0.995)+quantile(sh30[,3,q],0.995)+quantile(sh30[,4,q],0.995)+
                   quantile(sh30[,5,q],0.995)+quantile(sh30[,6,q],0.995)+quantile(sh30[,7,q],0.995)+quantile(sh30[,8,q],0.995)+
                   quantile(sh30[,9,q],0.995)+quantile(sh30[,10,q],0.995)+quantile(sh30[,11,q],0.995)+quantile(sh30[,12,q],0.995)+
                   quantile(sh30[,13,q],0.995)+quantile(sh30[,14,q],0.995)+quantile(sh30[,15,q],0.995)+quantile(sh30[,16,q],0.995)+
                   quantile(sh30[,17,q],0.995)+quantile(sh30[,18,q],0.995)+quantile(sh30[,19,q],0.995)+quantile(sh30[,20,q],0.995)+
                   quantile(sh30[,21,q],0.995)+quantile(sh30[,22,q],0.995)+quantile(sh30[,23,q],0.995)+quantile(sh30[,24,q],0.995)+
                   quantile(sh30[,25,q],0.995)+quantile(sh30[,26,q],0.995)+quantile(sh30[,27,q],0.995)+quantile(sh30[,28,q],0.995)+
                   quantile(sh30[,29,q],0.995)+quantile(sh30[,30,q],0.995))/30
}


v95=sqrt(30)*var500095-varh500095
v99=sqrt(30)*var500099-varh500099
v995=sqrt(30)*var5000995-varh5000995

a95=sqrt(30)*quantile(dax_log_xts6000,0.95)-12.13 #0.397
a99=sqrt(30)*quantile(dax_log_xts6000,0.99)-22.87 #0.747
a995=sqrt(30)*quantile(dax_log_xts6000,0.995)-25.55  #2.645

#Teststatistik.nicht signifikant
a95/sd(v95) #0.40 
a99/sd(v99) #0.37 
a995/sd(v995) #1.09
2*(1-pnorm(0.40)) #0.69
2*(1-pnorm(0.37)) #0.71
2*(1-pnorm(1.09)) #0.28

#h=10
sh10=array(0,c(599,10,5000))  #599=nrow, 10=ncol, 5000=n Wiederholung
for (q in 1:5000){
  for (k in (1:10)){
    for (i in (1:599)){
      sh10[i,k,q]=sum(resam[,q][(10*i-10+k):(10*i+k-1)])
    }
  }
}

for (q in 1:5000){
  varh5000951[q]= (quantile(sh10[,1,q],0.95)+quantile(sh10[,2,q],0.95)+quantile(sh10[,3,q],0.95)+quantile(sh10[,4,q],0.95)+
                    quantile(sh10[,5,q],0.95)+quantile(sh10[,6,q],0.95)+quantile(sh10[,7,q],0.95)+quantile(sh10[,8,q],0.95)+
                    quantile(sh10[,9,q],0.95)+quantile(sh10[,10,q],0.95))/10
  varh5000991[q]= (quantile(sh10[,1,q],0.99)+quantile(sh10[,2,q],0.99)+quantile(sh10[,3,q],0.99)+quantile(sh10[,4,q],0.99)+
                    quantile(sh10[,5,q],0.99)+quantile(sh10[,6,q],0.99)+quantile(sh10[,7,q],0.99)+quantile(sh10[,8,q],0.99)+
                    quantile(sh10[,9,q],0.99)+quantile(sh10[,10,q],0.99))/10
  varh50009951[q]= (quantile(sh10[,1,q],0.995)+quantile(sh10[,2,q],0.995)+quantile(sh10[,3,q],0.995)+quantile(sh10[,4,q],0.995)+
                     quantile(sh10[,5,q],0.995)+quantile(sh10[,6,q],0.995)+quantile(sh10[,7,q],0.995)+quantile(sh10[,8,q],0.995)+
                     quantile(sh10[,9,q],0.995)+quantile(sh10[,10,q],0.995))/10
}


v951=sqrt(10)*var500095-varh5000951
v991=sqrt(10)*var500099-varh5000991
v9951=sqrt(10)*var5000995-varh50009951

a951=sqrt(10)*quantile(dax_log_xts6000,0.95)-6.75 #0.4826
a991=sqrt(10)*quantile(dax_log_xts6000,0.99)-13.34 #0.2954
a9951=sqrt(10)*quantile(dax_log_xts6000,0.995)-16.16  #0.1182

#Teststatistik.nicht signifikant
a951/sd(v951) #1.53 
a991/sd(v991) #0.33 
a9951/sd(v9951) #0.10 
2*(1-pnorm(1.53)) #0.13
2*(1-pnorm(0.33)) #0.74
2*(1-pnorm(0.10)) #0.92

###Modifizierte Square-Root-of-Time
#h=10
mvar9510=sqrt(10*VR10)*quantile(dax_log_xts6000,0.95)  #6.83
mvar9910=sqrt(10*VR10)*quantile(dax_log_xts6000,0.99)  #12.88
mvar99510=sqrt(10*VR10)*quantile(dax_log_xts6000,0.995)  #15.38
##zeta(30)
(6.83/6.75-1)*100 #1.19
(12.88/13.34-1)*100 #-3.45
(15.38/16.16-1)*100 #-4.83

#subsample-based test
(6.83-6.75)/sd(v951)    #0.25
(12.88453-13.34)/sd(v991) #-0.50
(15.38185-16.16)/sd(v9951) #-0.65
2*(1-pnorm(0.25)) #0.80
2*pnorm(-0.50) #0.62
2*pnorm(-0.65)#0.52


#h=30
mvar9530=sqrt(30*VR10)*quantile(dax_log_xts6000,0.95)  #11.84
mvar9930=sqrt(30*VR10)*quantile(dax_log_xts6000,0.99)  #22.32
mvar99530=sqrt(30*VR10)*quantile(dax_log_xts6000,0.995)  #26.64
#zeta(30)
(11.84/12.13-1)*100 #-2.39
(22.32/22.87-1)*100 #-2.40
(26.64/25.55-1)*100 #4.27

#subsample-based test
(11.84-12.13)/sd(v951)    #-0.92
(22.32-22.87)/sd(v991) #-0.61
(26.64-25.55)/sd(v9951) #0.91
2*pnorm(-0.92)#0.38
2*pnorm(-0.61)#0.54
2*(1-pnorm(0.91)) #0.36

##nicht signifikant, bietet kleine Verbesserung. Die Ergebnisse decken sich mit Cheng et al.2011

###Tail-Index Based Skalierung (Vries 1998, Jon Danielsson 2011)


