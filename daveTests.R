library(nimble)
library(cutoff2)
rm(list=ls())

mu = c(4, 6)
sd = c(1, 4)
w  = c(0.5, 0.5) # nimble::rdirch(1, c(1,1)/2)
# Plot mixture
xRange = range(qnorm(p=c(0.001,0.001,0.999,0.999),mu,sd))
xRange[1] = floor(xRange[1])
xRange[2] = ceiling(xRange[2])
curve(w[1]*dnorm(x,mu[1],sd[1]) + w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111)
curve(w[1]*dnorm(x,mu[1],sd[1]), xRange[1], xRange[2], n=1111, add=TRUE, col="blue")
curve(w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111, add=TRUE, col="red")
print(w)

set.seed(11)
n = 100
n1 = rbinom(1,n,w[1]); n2 = n - n1
y = c(rnorm(n1,mu[1],sd[1]), rnorm(n2,mu[2],sd[2]))
miny = floor(min(y))
maxy = ceiling(max(y))
##
hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.2))
curve(w[1]*dnorm(x,mu[1],sd[1])+w[2]*dnorm(x,mu[2],sd[2]),miny,maxy,add=TRUE,col="black", lwd=3)
curve(w[1]*dnorm(x,mu[1],sd[1]),miny,maxy,add=TRUE,col="black", lwd=2, lty=2)
curve(w[2]*dnorm(x,mu[2],sd[2]),miny,maxy,add=TRUE,col="black", lwd=2, lty=2)

# Estimate parameters of finite mixture model:
(fit1 <- em(y,"normal","normal"))
#' # Adding the E-M estimated finite mixture model:
lines(fit1, col="tomato", lwd=2)
#' # Estimate a cutoff from the fitted mixture model (Choisy method):
(cut_off <- cutoff(fit1, whose="Titterington", nb=1000))
polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55), col=rgb(1, 0.39, 0.28,.2),border=NA)
abline(v=cut_off[-1],lty=2,col="tomato")
abline(v=cut_off[1],col="tomato")

#### Penalised fit
# Estimate parameters of finite mixture model:
## debug(calcPenalty)
(fit2 <- em(y,"normal","normal", penaltyScale=1E4))
#' # Adding the E-M estimated finite mixture model:
lines(fit2, col="red", lwd=2)
#' # Estimate a cutoff from the fitted mixture model (Choisy method):
(cut_off <- cutoff(fit2, whose="Titterington", nb=1000))
polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55), col=rgb(1,0,0,.2),border=NA)
abline(v=cut_off[-1],lty=2,col="red")
abline(v=cut_off[1],col="red")
