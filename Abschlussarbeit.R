# Packages 
#install.packages("extRemes")
#install.packages("xts")
library(extRemes) #Extremwerttheorie
library(xts)      #eXtensible Time Series

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


#Jaehrliche Maxima
jahrmax=apply.yearly(dax_log_xts,min) #die jaehliche groesste Verlust 
jahrmax=-jahrmax  #Votzeichen umkeheren, damit sie positiv sind
jahrmax
plot(jahrmax)



