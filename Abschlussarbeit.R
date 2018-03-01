# Packages 
#install.packages("extRemes")
#install.packages("xts")
library(extRemes) #Extremwerttheorie
library(xts)      #eXtensible Time Series
library(tseries)  #Zeitreihe
library(MASS)
library(rugarch)  #univariate GARCH
library(skewt)    #Skewed-t-Verteilung
library(fGarch)   #GARCH
library(forecast)
library(evir)   #auch Extremwerttheorie

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
plot(dax_log_xts)                                      #Plot von Log_Return, aber warum chinesisch?
length(dax_log_xts)#Laenge=6826
#plot(dax_log$date,dax_log$logreturn,type = "l")
hist(dax_log_xts,breaks=100)                       #Hist von Log-Return
stats::qqnorm(dax_log$logreturn);qqline(dax_log$logreturn)  #leptokurtisch

adf.test(dax_log_xts)   #ADF-Test: p=0.01, deshalb ist es stationaer.

dax_log_xts=-dax_log_xts  # Vorzeichen umkehren, Verlust statt Rendite. --Negative Returns
plot(dax_log_xts)

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
#bmm_monat_st=period.max(dax_log$logreturn,seq(from=1,to=6826,by=20))
#bmm_quartal_st=period.max(dax_log$logreturn,seq(from=1,to=6826,by=60))
#bmm_semester_st=period.max(dax_log$logreturn,seq(from=1,to=6826,by=120))
#bmm_jahr_st=period.max(dax_log$logreturn,seq(from=1,to=6826,by=240))



####Monat:
fit_monat_st <- fevd(as.vector(bmm_monat_st), method = "MLE", type="GEV") 
x_monat_st=sort(fit_monat_st$x)    ####von klein zu gross geordnet
z_monat_st=pgev(x_monat_st,xi=fit_monat_st$results$par[3],mu=fit_monat_st$results$par[1],sigma=fit_monat_st$results$par[2])
##Stephens 1977 S.2
#statistic_W_1=0
#for (i in (1:length(z_monat_st))){
#  statistic_W_1=statistic_W_1+(z_monat_st[i]-(2*i-1)/(2*length(z_monat_st)))^2
#}
#statistic_W=statistic_W_1+1/(12*length(z_monat_st))
#statistic_W_m=statistic_W*(1+0.2/sqrt(length(z_monat_st)))  ###0.33>0.124 #H0 abgelehnt.
###Goodness of Fit Shermann 1957
monat_order=sort(as.numeric(bmm_monat_st))
length(monat_order) #343
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
omega #statistik = 0.3892
mean_monat=(343/344)^(343+1)
var_monat=(2*exp(1)-5)/(exp(2)*343)
qnorm(0.95,mean=mean_monat,sd=sqrt(var_monat)) #Kritischer Wert= 0.3889 < 0.3892 H0 abgelehnt

###Quartal
fit_quartal_st <- fevd(as.vector(bmm_quartal_st), method = "MLE", type="GEV") 
x_quartal_st=sort(fit_quartal_st$x)    ####von klein zu gross geordnet
z_quartal_st=pgev(x_quartal_st,xi=fit_quartal_st$results$par[3],mu=fit_quartal_st$results$par[1],sigma=fit_quartal_st$results$par[2])
##Stephens 1977 S.2
#statistic_W_1=0
#for (i in (1:length(z_quartal_st))){
#  statistic_W_1=statistic_W_1+(z_quartal_st[i]-(2*i-1)/(2*length(z_quartal_st)))^2
#}
#statistic_W=statistic_W_1+1/(12*length(z_quartal_st))
#statistic_W_m=statistic_W*(1+0.2/sqrt(length(z_quartal_st)))  ###0.35>0.124 #H0 abgelehnt.
quartal_order=sort(as.numeric(bmm_quartal_st))
length(quartal_order) #115
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
omega #statistik = 0.4088
mean_quartal=(115/116)^(115+1)
var_quartal=(2*exp(1)-5)/(exp(2)*115)
qnorm(0.95,mean=mean_quartal,sd=sqrt(var_quartal)) #Kritischer Wert= 0.4035< 0.4088 H0 abgelehnt


#Semester
fit_semester_st <- fevd(as.vector(bmm_semester_st), method = "MLE", type="GEV") 
x_semester_st=sort(fit_semester_st$x)    ####von klein zu gross geordnet
z_semester_st=pgev(x_semester_st,xi=fit_semester_st$results$par[3],mu=fit_semester_st$results$par[1],sigma=fit_semester_st$results$par[2])
##Stephens 1977 S.2
#statistic_W_1=0
#for (i in (1:length(z_semester_st))){
#  statistic_W_1=statistic_W_1+(z_semester_st[i]-(2*i-1)/(2*length(z_semester_st)))^2
#}
#statistic_W=statistic_W_1+1/(12*length(z_semester_st))
#statistic_W_m=statistic_W*(1+0.2/sqrt(length(z_semester_st)))  ###0.086<0.124 #H0 nicht abgelehnt.
semester_order=sort(as.numeric(bmm_semester_st))
length(semester_order) #58
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
omega #statistik = 0.4124
mean_semester=(58/59)^(58+1)
var_semester=(2*exp(1)-5)/(exp(2)*58)
qnorm(0.95,mean=mean_semester,sd=sqrt(var_semester)) #Kritischer Wert= 0.4172>0.4124 H0 nicht ab


##Jahr
fit_jahr_st <- fevd(as.vector(bmm_jahr_st), method = "MLE", type="GEV") 
x_jahr_st=sort(fit_jahr_st$x)    ####von klein zu gross geordnet
z_jahr_st=pgev(x_jahr_st,xi=fit_jahr_st$results$par[3],mu=fit_jahr_st$results$par[1],sigma=fit_jahr_st$results$par[2])
##Stephens 1977 S.2
#statistic_W_1=0
#for (i in (1:length(z_jahr_st))){
#  statistic_W_1=statistic_W_1+(z_jahr_st[i]-(2*i-1)/(2*length(z_jahr_st)))^2
#}
#statistic_W=statistic_W_1+1/(12*length(z_jahr_st))
#statistic_W_m=statistic_W*(1+0.2/sqrt(length(z_jahr_st)))  ###0.04<0.124 #H0 nicht abgelehnt.
jahr_order=sort(as.numeric(bmm_jahr_st))
length(jahr_order) #30
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
omega #statistik = 0.36966
mean_jahr=(30/31)^(30+1)
var_jahr=(2*exp(1)-5)/(exp(2)*30)
qnorm(0.95,mean=mean_jahr,sd=sqrt(var_jahr)) #Kritischer Wert= 0.4348 > 0.36966 H0 nicht ab.


####semesterliche und jaehrliche bestehen bei dem Test. Deshalb wird n = 120 gewaehlt! (mehr Beobachtungen als die jaehlichen)  


# Moving Window (Groesse=2400) Zum Beispiel: das erste Fenster. Laenge=2400, damit durch 20,60,120,240 perfekt aufteilbar. 1200 zu klein fuer semesterlich.
ts_bm_1=dax_log_xts[1:2400]   
semestermax1=period.max(ts_bm_1,seq(from=120,to=2400,by=120))   #die groessete semesterliche Verlust
semestermax1
plot(semestermax1)
fit_semester1 <- fevd(as.vector(semestermax1), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
fit_semester1_evir=gev(dax_log_xts,block = 120)
plot(fit_semester1_evir)         #passt
fit_semester1$results$par  #Parameter extrahieren 

# VaRs und ESs berechnen. n=120 ist semesterlich
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
ES95_bmm=numeric(0)
ES99_bmm=numeric(0)
ES995_bmm=numeric(0)
VaR120 = function(x){mu-sigma/tau*(1-(-log((1-x)^(120)))^(-tau))} #tau = xi = shape Parameter
fit=numeric(0)
#returnlevel = function(x,mu,sigma,xi){mu-sigma/xi*(1-(-log((1-x)^120))^(-xi))}
for (i in (1:4426)){         #es gibt (6826-2400) Vorhersagen
  semestermax=period.max(dax_log_xts[i:(i+2399)],seq(from=120,to=2400,by=120))    #die groessete semesterliche Verlust
  fit <- fevd(as.vector(semestermax), method = "MLE", type="GEV")
  #fit = gev(dax_log_xts[i:(i+2399)],block=120)
  #VaR995_bmm[i]=returnlevel(0.005,mu=fit$par.ests[3],sigma=fit$par.ests[2],xi=fit$par.ests[1])#rlevel.gev(fit,2.212322)[2]
   # VaR99_bmm[i]=returnlevel(0.01,mu=fit$par.ests[3],sigma=fit$par.ests[2],xi=fit$par.ests[1])#rlevel.gev(fit,1.427308)[2]
  #  VaR95_bmm[i]=returnlevel(0.05,mu=fit$par.ests[3],sigma=fit$par.ests[2],xi=fit$par.ests[1])#rlevel.gev(fit,1.002127)[2]
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 2.212322)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^n=(1-1/k). Hier n = 120
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.427308)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.002127)[1]
   mu=as.numeric(fit$results$par[1])
    sigma=as.numeric(fit$results$par[2])
    tau=as.numeric(fit$results$par[3])
    ES995_bmm[i]=integrate(VaR120,lower=0,upper=0.005)$value*200
    ES99_bmm[i]=integrate(VaR120,lower=0,upper=0.01)$value*100
    ES95_bmm[i]=integrate(VaR120,lower=0,upper=0.05)$value*20
    }

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[2401:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[2401:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[2401:6826]) 
ES995_bmm_xts=xts(ES995_bmm,dax_log$date[2401:6826])
ES99_bmm_xts=xts(ES99_bmm,dax_log$date[2401:6826])
ES95_bmm_xts=xts(ES95_bmm,dax_log$date[2401:6826])

plot(dax_log_xts[2401:6826],main="VaR")  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")   
legend("bottomleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))


plot(dax_log_xts[2401:6826],main="ES")  
lines(ES995_bmm_xts,col="red")   
lines(ES99_bmm_xts,col="blue")    
lines(ES95_bmm_xts,col="green")  
legend("bottomleft",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_bmm_xts<dax_log_xts[2401:6826]) # 82 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[2401:6826]) #  152 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[2401:6826]) #  549 Ueberschreitungen


# Mit Theta (Extremaler Index) Embrechts 1998 Chap8 P.S19

#n=60,n=120. N_u=15,20,25,30,40,50,100,200

#n=60  #Mcneil 1998 S.13,14  (m, n muessen gross genug sein)
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

#n=120
n120max=period.max(dax_log_xts,seq(from=120,to=6826,by=120))
fit_120 <- fevd(as.vector(n120max), method = "MLE", type="GEV")  #MLE um Parameter zu schaetzen
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
K15=sum(n120max>N15);K20=sum(n120max>N20);K25=sum(n120max>N25);K30=sum(n120max>N30);K40=sum(n120max>N40);
K50=sum(n120max>N50);K100=sum(n120max>N100);K200=sum(n120max>N200)
K15;K20;K25;K30;K40;K50;K100;K200
theta=function(n,m,Ku,Nu){(log(1-Ku/m)/log(1-Nu/6826))/n}
#Theta berechnen
theta15=theta(n=120,m=57,Ku=K15,Nu=15)
theta20=theta(n=120,m=57,Ku=K20,Nu=20)
theta25=theta(n=120,m=57,Ku=K25,Nu=25)
theta30=theta(n=120,m=57,Ku=K30,Nu=30)
theta40=theta(n=120,m=57,Ku=K40,Nu=40)
theta50=theta(n=120,m=57,Ku=K50,Nu=50)
theta100=theta(n=120,m=57,Ku=K100,Nu=100)
theta200=theta(n=120,m=57,Ku=K200,Nu=200)
theta_dach=mean(c(theta15,theta20,theta25,theta30,theta40,theta50,theta100,theta200))#Mcneil 1998 S.13,14
theta_dach #0.42  



#Theta mir "evir" berechnen. n=60,120 ...

#n=120 Longin 2000: Block = Semester
n120max=period.max(dax_log_xts,seq(from=120,to=6826,by=120))
#sum(n120max>quantile(dax_log_xts,0.95))
sum(n120max>5) #Ueberschreiungen des Blocks = 13. Longin 2000: 5% als Threshold. Und Block = Semester
index120=exindex(dax_log_xts,block=120)
index120 #0.3190

#sum(n60max>5)  #17  #Longin 2000: 5% als Threshold. 
#index60=exindex(dax_log_xts,block=60) 
#index60 #0.398




#VaR
#n=120, Moving Window = 2400. Theta= 0.3190
VaR95_bmm=numeric(0)
VaR99_bmm=numeric(0)
VaR995_bmm=numeric(0)
ES95_bmm=numeric(0)
ES99_bmm=numeric(0)
ES995_bmm=numeric(0)
n120max=numeric(0)
fit=numeric(0)
VaR120 = function(x){mu-sigma/tau*(1-(-log((1-x)^(120*0.3190)))^(-tau))} #tau = xi = shape Parameter

for (i in (1:4426)){         #es gibt (6826-2400) Vorhersagen
  n120max=period.max(dax_log_xts[i:(i+2399)],seq(from=120,to=2400,by=120))    #die groessete quartalliche Verlust
  fit <- fevd(as.vector(n120max), method = "MLE", type="GEV")
  VaR995_bmm[i]=return.level(fit, conf = 0.05, return.period= 5.727568)[1] #Umrechnung zwischen r.p und Quantil, Siehe Longin2000, Mcneil1998. (1-p)^(n*theta)=(1-1/k). Hier n = 120, Theta=0.3190
  VaR99_bmm[i]=return.level(fit, conf = 0.05, return.period= 3.131228)[1]  
  VaR95_bmm[i]=return.level(fit, conf = 0.05, return.period= 1.163285)[1] 
  mu=as.numeric(fit$results$par[1])
  sigma=as.numeric(fit$results$par[2])
  tau=as.numeric(fit$results$par[3])
  ES995_bmm[i]=integrate(VaR120,lower=0,upper=0.005)$value*200
  ES99_bmm[i]=integrate(VaR120,lower=0,upper=0.01)$value*100
  ES95_bmm[i]=integrate(VaR120,lower=0,upper=0.05)$value*20
}

VaR995_bmm_xts=xts(VaR995_bmm,dax_log$date[2401:6826])
VaR99_bmm_xts=xts(VaR99_bmm,dax_log$date[2401:6826])
VaR95_bmm_xts=xts(VaR95_bmm,dax_log$date[2401:6826]) 
ES995_bmm_xts=xts(ES995_bmm,dax_log$date[2401:6826])
ES99_bmm_xts=xts(ES99_bmm,dax_log$date[2401:6826])
ES95_bmm_xts=xts(ES95_bmm,dax_log$date[2401:6826])

plot(dax_log_xts[2401:6826],main="VaR")  
lines(VaR995_bmm_xts,col="red")   
lines(VaR99_bmm_xts,col="blue")    
lines(VaR95_bmm_xts,col="green")  
legend("bottomleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))


plot(dax_log_xts[2401:6826],main="ES")  
lines(ES995_bmm_xts,col="red")   
lines(ES99_bmm_xts,col="blue")    
lines(ES95_bmm_xts,col="green")  
legend("bottomleft",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_bmm_xts<dax_log_xts[2401:6826]) #  23 Ueberschreitungen
sum(VaR99_bmm_xts<dax_log_xts[2401:6826]) #  57 Ueberschreitungen
sum(VaR95_bmm_xts<dax_log_xts[2401:6826]) #  231 Ueberschreitungen



###############POT################
# Mean Residual Life Plot: (Mean Excess)
mrlplot(dax_log_xts, main="Mean Residual Life Plot")    #u ist vielleicht in (0,4), aber nicht informative
meplot(dax_log_xts,xlim=c(0,5),ylim=c(1,1.5),type="l")  #u ist vielleicht 3.5. Nach 3.5 ist linear. Ist (6826-100)/6826=0.9854Quantil


#Hill-Schaetzer (Hill, 1975) (Mecneil 2000) tau^dach=1/k*Sigma^k_(j=1)(log(zj)-log(z(k+1))). 

#####Aber bei Hill-Schaetzer: shape-Parameter muss >0! (Mcneil 2000 Seite.17) Kann noch als Instrument zur Threshold-Wahl?
n=6826
evir::hill(dax_log_xts,xlim=c(15,300))  #Hill Plot.Ab 100 ist es linear  k=ungefaehr 100, y-Achse = ungefaehr 3.2.
quantile(dax_log_xts,(6826-100)/6826)  # Threshold wird als 3.50 gewaehlt.


#dax_log_order=sort(-dax_log$logreturn,decreasing = TRUE) 
#dax_log_order[100]#u = 3.50

taudach=numeric(0)
for (i in (15:600)){
  taudach[i]=1/i*sum(log(dax_log_order[1:i])-log(dax_log_order[i]))
}
plot(taudach,type="l")    #identisch zur evir::hill. Ab ungefaehr 100 ise es linear

# mit unterschiedlichen Grenzwerten
threshrange.plot(dax_log_xts, r = c(0, 5), nint = 16)
# ismev Implementation ist schneller:
ismev::gpd.fitrange(dax_log_xts, umin=0, umax=5, nint = 16) 


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
  #Mcneil,2000 S.17 the choice of k in Moving Window (10%)
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

plot(dax_log_xts[2401:6826],main="VaR")  
lines(VaR995_pot_xts,col="red")   
lines(VaR99_pot_xts,col="blue")    
lines(VaR95_pot_xts,col="green") 
legend("bottomleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))

plot(dax_log_xts[2401:6826],main="ES")  
lines(ES995_pot_xts,col="red")   
lines(ES99_pot_xts,col="blue")    
lines(ES95_pot_xts,col="green") 
legend("bottomleft",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))

sum(VaR995_pot_xts<dax_log_xts[2401:6826]) #  26 Ueberschreitungen
sum(VaR99_pot_xts<dax_log_xts[2401:6826]) #  61 Ueberschreitungen
sum(VaR95_pot_xts<dax_log_xts[2401:6826]) #  222 Ueberschreitungen



######Danielsson2001 Threshold mit Hilfe von Subsample-Bootstrap
#install.packages("tea")
library(tea) # Package zum Danielssons Bootstrap
danielsson(dax_log$logreturn,B=100) #Threshold= 4.20. aber das einzelne Verfahren kostet mehr als 20 Minuten

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
garchfit3@fit$par
garchfit3@residuals
garchfit3@sigma.t    #sd
garchfit3@h.t       #Var  
garchfit3@fitted    

zt=(dax_log_xts-garchfit3@fitted)/garchfit3@sigma.t   #standardisierte Residuen
plot(zt)
plot(abs(dax_log_xts))
plot(garchfit3@sigma.t,type="l")

par(mfrow=c(2,2))
acf(dax_log_xts);acf(abs(dax_log_xts))  #nicht i.i.d
acf(zt);acf(abs(zt))     #keinen ARCH-Effekt

Box.test(dax_log_xts,lag=10,type="Ljung-Box") #p=0.00 H0: unabhaengig verteilt
Box.test(zt,lag=10,type="Ljung-Box")  #p=0.53

stats::qqnorm(zt);qqline(dax_log$logreturn)

# Mean Residual Life Plot: (Mean Excess)
mrlplot(zt, main="Mean Residual Life Plot")    
meplot(zt,type="l")  
meplot(zt,type="l")  #Threshold = ungefaehr 2.5. Ab 2.5 ist es linear


#Hill-Plot
evir::hill(zt,xlim=c(15,500))  #Hill Plot  k=ungefaehr 70, y = ungefaehr 5.
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

#Moving Windows mit Laenge 2400. (Aber hier, die Parameter von GARCH veraendern sich nicht mit dem Moving Window)
VaR95_pot_garch=numeric(0)
VaR99_pot_garch=numeric(0)
VaR995_pot_garch=numeric(0)
ES95_pot_garch=numeric(0)
ES99_pot_garch=numeric(0)
ES995_pot_garch=numeric(0)
for (i in (1:4426)){         #es gibt (6826-2400) Vorhersagen. Threshold ist wie immer 90% Quantil
  gpdpotgarch=fevd(as.vector(zt[i:(2399+i)]), method = "MLE", type="GP", threshold=quantile(zt[i:(2399+i)],0.9))
  VaR995_pot_garch[i]=quantile(zt[i:(2399+i)],0.9)+gpdpotgarch$results$par[1]/gpdpotgarch$results$par[2]*((2400*0.005/240)^(-gpdpotgarch$results$par[2])-1)
  VaR99_pot_garch[i]=quantile(zt[i:(2399+i)],0.9)+gpdpotgarch$results$par[1]/gpdpotgarch$results$par[2]*((2400*0.01/240)^(-gpdpotgarch$results$par[2])-1)
  VaR95_pot_garch[i]=quantile(zt[i:(2399+i)],0.9)+gpdpotgarch$results$par[1]/gpdpotgarch$results$par[2]*((2400*0.05/240)^(-gpdpotgarch$results$par[2])-1)
}  

VaR995_pot_xts_garch=xts(VaR995_pot_garch*garchfit3@sigma.t[2401:6826]+garchfit3@fitted[2401:6826],dax_log$date[2401:6826])#Mcneil S.6
VaR99_pot_xts_garch=xts(VaR99_pot_garch*garchfit3@sigma.t[2401:6826]+garchfit3@fitted[2401:6826],dax_log$date[2401:6826])
VaR95_pot_xts_garch=xts(VaR95_pot_garch*garchfit3@sigma.t[2401:6826]+garchfit3@fitted[2401:6826],dax_log$date[2401:6826]) 

plot(dax_log_xts[2401:6826])  
lines(VaR995_pot_xts_garch,col="red")   
lines(VaR99_pot_xts_garch,col="blue")    
lines(VaR95_pot_xts_garch,col="green")   
sum(VaR995_pot_xts_garch<dax_log_xts[2401:6826]) #  17 Ueberschreitungen
sum(VaR99_pot_xts_garch<dax_log_xts[2401:6826]) #  34 Ueberschreitungen
sum(VaR95_pot_xts_garch<dax_log_xts[2401:6826]) #  234 Ueberschreitungen

#ES berechnen  Mcneil 2000 S.23 Gleichung-(17) 
ES995_pot_xts_garch=xts(garchfit3@fitted[2401:6826]+garchfit3@sigma.t[2401:6826]*VaR995_pot_garch*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(zt[i:(2399+i)],0.9))/(VaR995_pot_garch-VaR995_pot_garch*gpdpotgarch$results$par[2])),dax_log$date[2401:6826])
ES99_pot_xts_garch=xts(garchfit3@fitted[2401:6826]+garchfit3@sigma.t[2401:6826]*VaR99_pot_garch*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(zt[i:(2399+i)],0.9))/(VaR99_pot_garch-VaR995_pot_garch*gpdpotgarch$results$par[2])),dax_log$date[2401:6826])
ES95_pot_xts_garch=xts(garchfit3@fitted[2401:6826]+garchfit3@sigma.t[2401:6826]*VaR95_pot_garch*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(zt[i:(2399+i)],0.9))/(VaR95_pot_garch-VaR995_pot_garch*gpdpotgarch$results$par[2])),dax_log$date[2401:6826])


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
  #ES: Mcneil S.23
  ES995_pot_garch[i]=predict(garchfitm)[1,1]+predict(garchfitm)[1,3]*VaR995_pot_z*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(ztm,0.9))/(VaR995_pot_z-VaR995_pot_z*gpdpotgarch$results$par[2]))
  ES99_pot_garch[i]=predict(garchfitm)[1,1]+predict(garchfitm)[1,3]*VaR99_pot_z*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(ztm,0.9))/(VaR99_pot_z-VaR995_pot_z*gpdpotgarch$results$par[2]))
  ES95_pot_garch[i]=predict(garchfitm)[1,1]+predict(garchfitm)[1,3]*VaR95_pot_z*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(ztm,0.9))/(VaR95_pot_z-VaR995_pot_z*gpdpotgarch$results$par[2]))
}  
VaR995_pot_xts_garch=xts(VaR995_pot_garch,dax_log$date[2401:6826])
VaR99_pot_xts_garch=xts(VaR99_pot_garch,dax_log$date[2401:6826])
VaR95_pot_xts_garch=xts(VaR95_pot_garch,dax_log$date[2401:6826])


plot(dax_log_xts[2401:6826],main="VaR")  
lines(VaR995_pot_xts_garch,col="red")   
lines(VaR99_pot_xts_garch,col="blue")    
lines(VaR95_pot_xts_garch,col="green")   
legend("bottomleft",inset=0.005,c("VaR0.995","VaR0.99","VaR0.95"),col=c("red","blue","green"),lty=c(1,1,1))
sum(VaR995_pot_xts_garch<zt[2401:6826]) #  28 Ueberschreitungen
sum(VaR99_pot_xts_garch<zt[2401:6826]) #  43 Ueberschreitungen
sum(VaR95_pot_xts_garch<zt[2401:6826]) #  199 Ueberschreitungen


#ES berechnen  Mcneil 2000 S.23 Gleichung-(17) 

ES995_pot_xts_garch=xts(garchfit3@fitted[2401:6826]+garchfit3@sigma.t[2401:6826]*VaR995_pot_garch*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(zt[i:(2399+i)],0.9))/(VaR995_pot_garch-VaR995_pot_garch*gpdpotgarch$results$par[2])),dax_log$date[2401:6826])
ES99_pot_xts_garch=xts(garchfit3@fitted[2401:6826]+garchfit3@sigma.t[2401:6826]*VaR99_pot_garch*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(zt[i:(2399+i)],0.9))/(VaR99_pot_garch-VaR995_pot_garch*gpdpotgarch$results$par[2])),dax_log$date[2401:6826])
ES95_pot_xts_garch=xts(garchfit3@fitted[2401:6826]+garchfit3@sigma.t[2401:6826]*VaR95_pot_garch*(1/(1-gpdpotgarch$results$par[2])+(gpdpotgarch$results$par[1]-gpdpotgarch$results$par[2]*quantile(zt[i:(2399+i)],0.9))/(VaR95_pot_garch-VaR995_pot_garch*gpdpotgarch$results$par[2])),dax_log$date[2401:6826])
plot(dax_log_xts[2401:6826],main="ES")  
lines(ES995_pot_xts_garch,col="red")   
lines(ES99_pot_xts_garch,col="blue")    
lines(ES95_pot_xts_garch,col="green") 
legend("bottomleft",inset=0.005,c("ES0.995","ES0.99","ES0.95"),col=c("red","blue","green"),lty=c(1,1,1))




