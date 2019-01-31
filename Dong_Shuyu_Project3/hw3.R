library(ggplot2)

Greeks=read.csv("Shuyu_greeks.csv",header=F)
Greeks=Greeks[,-8]
colnames(Greeks)=c("S0","Call_price","Delta","Theta","Vega","Gamma","Rho")



par(mfrow=c(2,3))
plot(y=Greeks$Delta,x=Greeks$S0,main="Delta",type = "l")
plot(y=Greeks$Theta,x=Greeks$S0,main="Theta",type = "l")
plot(y=Greeks$Vega,x=Greeks$S0,main="Vega",type="l")
plot(y=Greeks$Gamma,x=Greeks$S0,main="Gamma",type="l")
plot(y=Greeks$Rho,x=Greeks$S0,main="Rho",type="l")

halton=read.csv("Shuyu_Haltons.csv",header = F)

par(mfrow=c(2,2))

plot(halton$V1,halton$V2,pch=20,main="LGM Uniform")
plot(halton$V3,halton$V4,pch=20,main="Halton bases(2,7)")
plot(halton$V5,halton$V6,pch=20,main="Halton bases(2,4)")

