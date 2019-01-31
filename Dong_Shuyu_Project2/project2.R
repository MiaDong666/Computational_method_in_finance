library(xts)
library(ggplot2)

stock1=read.csv("Shuyu_Stock1.csv",header = F)
stockpath1=read.csv("Shuyu_Stockpath1.csv",header=F)

stock1=stock1[,-2]

plot(stock1,main="E[Sn] sigma=0.18",pch=20)


stockpath1=stockpath1[,-7]
stockpath1$X=c(0:1000)
stockpath1$X=stockpath1$X/100
plot(stockpath1$V1,x=stockpath1$X,type="l",col="blue",ylim=c(0,450),main="stockpath & E[Sn] sigma=0.18")
points(stockpath1$V2,x=stockpath1$X,type="l",col="red")
points(stockpath1$V3,x=stockpath1$X,type="l",col="violet")
points(stockpath1$V4,x=stockpath1$X,type="l",col="green")
points(stockpath1$V5,x=stockpath1$X,type="l",col="orange")
points(stockpath1$V6,x=stockpath1$X,type = "l",col="cyan")
points(stock1,col="black",pch=20)
points(stock1,col="black",type="l")

stock2=read.csv("Shuyu_Stock2.csv",header = F)
stockpath2=read.csv("Shuyu_Stockpath2.csv",header=F)

stock2=stock2[,-2]

plot(stock2,main="E[Sn] sigma=0.35",pch=20)


stockpath2=stockpath2[,-7]
stockpath2$X=c(0:1000)
stockpath2$X=stockpath2$X/100
plot(stockpath2$V1,x=stockpath2$X,type="l",col="blue",ylim=c(0,600),main="stockpath & E[Sn] sigma=0.35")
points(stockpath2$V2,x=stockpath2$X,type="l",col="red")
points(stockpath2$V3,x=stockpath2$X,type="l",col="violet")
points(stockpath2$V4,x=stockpath2$X,type="l",col="green")
points(stockpath2$V5,x=stockpath2$X,type="l",col="orange")
points(stockpath2$V6,x=stockpath2$X,type = "l",col="cyan")
points(stock2,col="black",pch=20)
points(stock2,col="black",type="l")  



