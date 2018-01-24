daten=read.csv("DAX.csv") #Quelle:finance.yahoo.com
daten_omit=daten[!daten$Open=="null",] #Fehlende Werte weg
dim(daten_omit) #noch 6827 Beobachtungen
dax=daten_omit[,1:2]#Nur Datum und Schlusskurs bleiben
dax$Close=as.numeric(as.character(dax$Close))
plot(dax)

n=length(dax$Close)
logreturn <- 100*log(dax$Close[-1]/dax$Close[-n])#log Return
dax_log=data.frame(dax$Date[-1],logreturn)
names(dax_log)[names(dax_log)=="dax.Date..1."] <- "date"  #Umbennen
head(dax_log)
plot(as.Date(dax_log$date),dax_log$logreturn,type = "l")
hist(dax_log$logreturn,breaks=100)
