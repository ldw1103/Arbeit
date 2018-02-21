library(ggplot2)
library(MASS)
library(grid)#Lokation von ggplots2
library(fGarch)#Skewed-t-Verteilung
library(skewt)#Skewed-t-Verteilung

#Seite 12. Studentsche t-Verteilung und NV-Verteilung
tfitt=fitdistr(-dax_log$logreturn,"t")
tfitt
m=tfitt$estimate[1]
s=tfitt$estimate[2]
df=tfitt$estimate[3]
pi=3.1415926
ft=function(x,df,m,s){(gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5))}

t_nv=ggplot(data=dax_log,aes(-logreturn))+ labs(x = "Verlust",y="Dichte")+
geom_histogram(aes(y=..density..),fill="blue",bins=100)+
labs(title="Histogramm des Log-Verlusts")+
#geom_line(aes(y = ..density..), stat = 'density',lwd=1.2)+
  stat_function(fun = dnorm, n = 1000, args = list(mean = mean(-dax_log$logreturn), sd = sd(-dax_log$logreturn)),lwd=1.1) +
  stat_function(fun = function(x) (gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5)) ,lwd=1.1,col="red") +
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),)
t_nv


t_nv1=ggplot(data=dax_log,aes(-logreturn))+ labs(x = "Verlust",y="Dichte")+coord_cartesian(xlim=c(3, 8),ylim=c(0, 0.05))+
  geom_histogram(aes(y=..density..),fill="blue",bins=200)+
  labs(title="Histogramm des Log-Verlusts")+
  #geom_line(aes(y = ..density..), stat = 'density',lwd=1.1,col="green")+
  #stat_function(fun = dnorm, n = 1000, args = list(mean = mean(-dax_log$logreturn), sd = sd(-dax_log$logreturn)),lwd=1.1) +
  stat_function(fun = function(x) (gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5)) ,lwd=1.1,col="red") +
  stat_function(fun = dskt, n = 1000, args = list(mean = mean(-dax_log$logreturn), sd = sd(-dax_log$logreturn)),lwd=1.1)+
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

skt=ggplot(data=dax_log,aes(-logreturn))+ labs(x = "Verlust",y="Dichte") +#coord_cartesian(xlim=c(3, 8),ylim=c(0, 0.05))+
  geom_histogram(aes(y=..density..),fill="blue",bins=200)+
  labs(title="Histogramm des Log-Verlusts")+
  #geom_line(aes(y = ..density..), stat = 'density',lwd=1.1,col="green")+
  stat_function(fun = dnorm, n = 1000, args = list(mean = mean(-dax_log$logreturn), sd = sd(-dax_log$logreturn)),lwd=1.1) +
  stat_function(fun = function(x) (gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5)) ,lwd=1.1,col="red") +
  stat_function(fun = dskt, n = 1000, args = list(df=sktfit$estimate[3],gamma=sktfit$estimate[4]),lwd=1.1,col="green") +
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),)
skt

skt1=ggplot(data=dax_log,aes(-logreturn))+ labs(x = "Verlust",y="Dichte") +coord_cartesian(xlim=c(3, 8),ylim=c(0, 0.05))+
  geom_histogram(aes(y=..density..),fill="blue",bins=200)+
  labs(title="Histogramm des Log-Verlusts")+
  #geom_line(aes(y = ..density..), stat = 'density',lwd=1.1,col="green")+
  stat_function(fun = dnorm, n = 1000, args = list(mean = mean(-dax_log$logreturn), sd = sd(-dax_log$logreturn)),lwd=1.1) +
  stat_function(fun = function(x) (gamma((df+1)/2))/((gamma(df/2)*sqrt(pi*df)*s))*(1+(1/df)*((x-m)/s)^2)^(-(0.5*df+0.5)) ,lwd=1.1,col="red") +
  stat_function(fun = dskt, n = 1000, args = list(df=sktfit$estimate[3],gamma=sktfit$estimate[4]),lwd=1.1,col="green") +
  theme(axis.text=element_text(size=20), axis.title.x = element_text(size=25),
        axis.title.y=element_text(size=25,face="bold"),title =element_text(size=20, face='bold'),)
skt1

print(skt,vp=vp1)
print(skt1,vp=vp2)


dsstd()

