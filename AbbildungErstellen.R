library(ggplot2)
library(MASS)
library(grid)#Lokation von ggplots2
library(fGarch)#Skewed-t-Verteilung
library(skewt)#Skewed-t-Verteilung
library(reshape2)

#Seite 12. Studentsche t-Verteilung und NV-Verteilung
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

##location
vp1= viewport(width=0.5,height=1,x = 0, y = 1, just = c("left", "top"))
vp2= viewport(width=0.5,height=1,x = 0.75, y = 1, just = c("center", "top"))

print(t_nv,vp=vp1)
print(t_nv1,vp=vp2)


#Abbildung zur skewed t Verteilung

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

#Seite 18 3 EVT Dichte

evtd=ggplot(data.frame(x=c(-4, 4)), aes(x)) +labs(x = "x",y="Dichte")  +
  stat_function(fun=devd,args = list(loc = 0, scale = 1, shape = 0, threshold = 0, log = F,type="GEV"),lwd=1.1)+
  stat_function(fun=devd,args = list(loc = 1, scale = 1, shape = 1, threshold = 0, log = F,type="GEV"),lwd=1.1,col="red")+
  stat_function(fun=devd,args = list(loc = -1, scale = 1, shape = -1, threshold = 0, log = F,type="GEV"),lwd=1.1,col="green")+
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),)
evtd


###evt mit legend

x  <- seq(-5, 5, 0.01) 
Gumbel <- devd(x,loc = 0, scale = 1, shape = 0, threshold = 0, log = F,type="GEV")
Frechet <- devd(x,loc = 1, scale = 1, shape = 1, threshold = 0, log = F,type="GEV")
Weibull <- devd(x,loc = -1, scale = 1, shape = -1, threshold = 0, log = F,type="GEV")
df <- data.frame(x,Gumbel,Frechet,Weibull)

#ggplot(df, aes(x)) +    
#  geom_line(aes(y=Gumbel),lwd=1.1) +
#  geom_line(aes(y=Frechet), colour="blue",lwd=1.1) +  
#  geom_line(aes(y=Weibull), colour="red",lwd=1.1)

gg <- melt(df,id="x")
ggplot(gg, aes(x=x, y=value, color=variable))+
  geom_line(lwd=1.1)+labs(x = "x",y="Dichte")+
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),legend.text=element_text(size=17),legend.title=element_blank())
  