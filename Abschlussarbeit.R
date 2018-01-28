# Packages 
#install.packages("extRemes")
#install.packages("xts")
library(extRemes) #Extremwerttheorie
library(xts)      #eXtensible Time Series
library(tseries)
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
stats::qqnorm(dax_log$logreturn);qqline(dax_log$logreturn)  #leptokurtotic

adf.test(dax_log_xts)   #ADF-Test: p=0.01, deshalb ist es stationaer.

#Jaehrliche Maxima
jahrmax=apply.yearly(dax_log_xts,min) #die jaehliche groesste Verlust 
jahrmax=-jahrmax  #Vorzeichen umkeheren, damit sie positiv sind
jahrmax
plot(jahrmax)

# maximum-likelihood fitting of the GEV distribution
fit_mle <- fevd(as.vector(jahrmax), method = "MLE", type="GEV")
# diagnostic plots
plot(fit_mle)
fit_mle$results$par     #Paramter. Location=3.55 (mu), Scale=1.39 (sigma), Shape=0.195 (xi) 
#xi>0, somit ist es Frechet-Verteilung

# return levels:
rl_mle <- return.level(fit_mle, conf = 0.05, return.period= c(2,5,10,20,50,100))

#VaR berechnen, µ«ÊÇ²»¹»
VaR_0.95=3.5501838-1.3903295/0.1951717*(1-(-log(0.95))^(-0.1951717))
VaR_0.95    # gleich 20-year level

# mit Theta (extremer Index) Embrechtschap7 P.289





# fitting of GEV distribution based on L-moments estimation
fit_lmom <- fevd(as.vector(jahrmax), method = "Lmoments", type="GEV")
# diagnostic plots
plot(fit_lmom)
# return levels:
rl_lmom <- return.level(fit_lmom, conf = 0.05, return.period= c(2,5,10,20,50,100))

# return level plots
par(mfcol=c(1,2))
# return level plot w/ MLE
plot(fit_mle, type="rl",
     main="Return Level Plot for Baernkopf w/ MLE",
     ylim=c(0,200), pch=16)
loc <- as.numeric(return.level(fit_mle, conf = 0.05,return.period=100))
segments(100, 0, 100, loc, col= 'midnightblue',lty=6)
segments(0.01,loc,100, loc, col='midnightblue', lty=6)

# return level plot w/ LMOM
plot(fit_lmom, type="rl",
     main="Return Level Plot for Baernkopf w/ L-Moments",
     ylim=c(0,200))
loc <- as.numeric(return.level(fit_lmom, conf = 0.05,return.period=100))
segments(100, 0, 100, loc, col= 'midnightblue',lty=6)
segments(0.01,loc,100, loc, col='midnightblue', lty=6)

# comparison of return levels
results <- t(data.frame(mle=as.numeric(rl_mle),
                        lmom=as.numeric(rl_lmom)))
colnames(results) <- c(2,5,10,20,50,100)
round(results,1)


