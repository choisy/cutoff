## library(nimble)
library(cutoff2)
rm(list=ls())

set.seed(4)
par(mfrow=c(2,2))
mu = c(4, 6)
sd = c(1, 4)
w  = c(0.5, 0.5)

#######################
## Plot 1: the model ##
#######################
xRange = range(qnorm(p=c(0.001,0.001,0.999,0.999),mu,sd))
xRange[1] = floor(xRange[1])
xRange[2] = ceiling(xRange[2])
curve(w[1]*dnorm(x,mu[1],sd[1]) + w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111,
      lwd=2, ylab="Density" ,xlab="MFI")
curve(w[1]*dnorm(x,mu[1],sd[1]), xRange[1], xRange[2], n=1111, lty=2, lwd=2 ,add=TRUE, col="blue")
curve(w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111, lty=2, lwd=2, add=TRUE, col="red")
title("Biological nonsense model")
legend(legend=c("model", "pos","neg"),"topright", lwd=c(2,2,2), col=c("black","red","blue"),
       bty="n" ,lty=c(1,2,2))
print(w)

######################
## Plot 2: the data ##
######################
n = 100
n1 = rbinom(1,n,w[1]); n2 = n - n1
y = c(rnorm(n1,mu[1],sd[1]), rnorm(n2,mu[2],sd[2]))
miny = floor(min(y))
maxy = ceiling(max(y))
hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.5), xlab="MFI", main="Simulated data")
curve(w[1]*dnorm(x,mu[1],sd[1])+w[2]*dnorm(x,mu[2],sd[2]),miny,maxy,add=TRUE,col="black", lwd=3)
curve(w[1]*dnorm(x,mu[1],sd[1]),miny,maxy,add=TRUE,col="black", lwd=2, lty=2)
curve(w[2]*dnorm(x,mu[2],sd[2]),miny,maxy,add=TRUE,col="black", lwd=2, lty=2)

#######################
## Unconstrained Fit ##
#######################
# Estimate parameters of finite mixture model:
(fit1 <- em(y,"normal","normal"))
hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.5), xlab="MFI", main="Unconstrained fit")
# Add the EM estimated finite mixture model:
lines(fit1, col="tomato", lwd=2)
# Estimate a cutoff from the fitted mixture model
(cut_off <- cutoff(fit1, whose="Titterington", nb=1000))
polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55), col=rgb(1, 0.39, 0.28,.2),border=NA)
abline(v=cut_off[-1],lty=2,col="tomato")
abline(v=cut_off[1],col="tomato")

###################
## Penalised fit ##
###################
# Estimate parameters of finite mixture model:
(fit2 <- em(y,"normal","normal", penaltyScale=1E4))
# Replot data
hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.5), xlab="MFI", main="Penalised fit")
# Add the penalised-EM estimated finite mixture model:
lines(fit2, col="red", lwd=2)
# Estimate a cutoff from the fitted mixture model
(cut_off <- cutoff(fit2, whose="Titterington", nb=1000))
polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55), col=rgb(1,0,0,.2),border=NA)
abline(v=cut_off[-1],lty=2,col="red")
abline(v=cut_off[1],col="red")
