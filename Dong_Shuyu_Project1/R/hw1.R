library(ggplot2)

#uniform
Uniform=read.csv("Shuyu_Uniform.csv",header = F)
Uniform=Uniform[,-3]
colnames(Uniform)=c("LGM","built-in Uniform")
#X
X=read.csv("Shuyu_X.csv",header = F)
X=X[,-2]
#binomial
binom=read.csv("Shuyu_binomial.csv",header = F)
binom=binom[,-2]
#expotential
expo=read.csv("Shuyu_exp.csv",header = F)
expo=expo[,-2]
#normal
normal=read.csv("Shuyu_Normal.csv",header = F)
normal=normal[,-3]
colnames(normal)=c("Box-muller","Polar-Marsaglia")


par(mfrow=c(1,2))
hist(Uniform$LGM,breaks=80,col="violet",main = "LGM")
hist(Uniform$`built-in Uniform`,breaks=80,col="blue",main="built-in function")

par(mfrow=c(1,1))
hist(X,breaks=80,col="violet")

hist(binom,breaks = 80,col="violet")

hist(expo,breaks=80,col="violet")

par(mfrow=c(1,2))
hist(normal$`Box-muller`,breaks=80,col="violet",main="Box-Muller")
hist(normal$`Polar-Marsaglia`,breaks=80,col="blue",main="Polar Marsaglia")





