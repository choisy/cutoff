library(nimble)
rm(list=ls())

mu = c(6, 4)
sd = c(1, 2)
w  = c(0.5, 0.5) # nimble::rdirch(1, c(1,1)/2)
# Plot mixture
xRange = range(qnorm(p=c(0.001,0.001,0.999,0.999),mu,sd))
xRange[1] = floor(xRange[1])
xRange[2] = ceiling(xRange[2])
curve(w[1]*dnorm(x,mu[1],sd[1]) + w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111)
curve(w[1]*dnorm(x,mu[1],sd[1]), xRange[1], xRange[2], n=1111, add=TRUE, col="blue")
curve(w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111, add=TRUE, col="red")
print(w)

dHash

cutoff(object,t=1e-64,nb=10,distr=2,type1=.05,level=.95,whose="Choisy") {
cutoff(mu1=mu[1],sigma1=sd[1],mu2=mu[2],sigma2=sd[2],w[1],D1="normal",D2="normal",distr=2,type1=.05,whose="Titterington")
c
## curve(w[1]*pnorm(x,mu[1],sd[1]) + w[2]*pnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111)
## curve(w[1]*pnorm(x,mu[1],sd[1]), xRange[1], xRange[2], n=1111, add=TRUE, col="blue")
## curve(w[2]*pnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111, add=TRUE, col="red")
## nimPrint(w)


## qTest = 1E-11
## (testQuantNeg = qnorm(p=qTest,mu[1],sd[1]))
## (testQuantPos = qnorm(p=1-qTest,mu[2],sd[2]))
## # Test in left hand tail
## w[1]*pnorm(testQuantNeg,mu[1],sd[1]) > w[2]*pnorm(testQuantNeg,mu[2],sd[2]) # TRUE indicates d(NegativeCurve) > d(PositiveCurve)
## # Test in right hand tail
## w[1]*pnorm(testQuantPos,mu[1],sd[1]) < w[2]*pnorm(testQuantPos,mu[2],sd[2]) # TRUE indicates d(PositiveCurve) > d(NegativeCurve)


penaltyScale = 10000
#### Penalise LHS
X = seq(qnorm(p=1E-11,mu[1],sd[1]), qnorm(p=1-1E-11,mu[1],sd[1]), l=1111)
# 1 - sum(dnorm(X,mu[1],sd[1])*delta) ## Almost zero, so sum is good approximation to one
iNonsenseLeft = w[1]*dnorm(X,mu[1],sd[1]) < w[2]*dnorm(X,mu[2],sd[2])
iNonsenseLeft = cumprod(iNonsenseLeft)==1
(ii = sum(iNonsenseLeft))
penaltyLeft = 0
if (ii > 0) {
  abline(v=X[ii], col="purple")
  (areaNonsenseLeft = (w[2]*(pnorm(X[ii],mu[2],sd[2]))) - (w[1]*(pnorm(X[ii],mu[1],sd[1]))))
  (penaltyLeft = log(1-areaNonsenseLeft))
}

#### Penalise RHS
X = seq(qnorm(p=1E-11,mu[2],sd[2]), qnorm(p=1-1E-11,mu[2],sd[2]), l=1111)
iNonsenseRight = w[1]*dnorm(X,mu[1],sd[1]) > w[2]*dnorm(X,mu[2],sd[2])
iNonsenseRight = rev(cumprod(rev(iNonsenseRight))==1)
ii = sum(!iNonsenseRight)
penaltyRight = 0
if (ii < length(X)) {
  abline(v=X[ii], col="orange")
  (areaNonsenseRight = (w[1]*(1-pnorm(X[ii],mu[1],sd[1]))) - (w[2]*(1-pnorm(X[ii],mu[2],sd[2])) ))
  (penaltyRight =  log(1-areaNonsenseRight))
}

#### Total Penalty
(penaltyScale * penaltyLeft)
(penaltyScale * penaltyRight)
(penaltyTotal = penaltyScale * (penaltyLeft + penaltyRight))

## ADD THIS TO MARC CHOISY'S LOG LIKELIHOOD

###########################
## Examples from package ##
###########################
sessionInfo()
library(mc2d)
library(MASS)
library(tree)

detach("package:cutoff", unload=TRUE)
source("~/R/package_workshop/cutoff/R/cutoff.R")
source("~/R/package_workshop/cutoff/R/em.R")
source("~/R/package_workshop/cutoff/R/estim_ci.R")
data(measles, package="cutoff")

# Measles IgG concentration data:
length(measles)
range(measles)
# Plotting the data:
hist(measles,100,F,xlab="concentration",ylab="density",ylim=c(0,.55), main=NULL,col="grey")
# The kernel density:
lines(density(measles),lwd=1.5,col="blue")
# Estimating the parameters of the finite mixture model:
(measles_out <- em(measles,"normal","normal"))
# The confidence interval of the parameter estimates:
confint(measles_out,t=1e-64,nb=100,level=.95)
# Adding the E-M estimated finite mixture model:
lines(measles_out,lwd=1.5,col="red")
# The legend:
legend("topleft",leg=c("non-parametric","E-M"),col=c("blue","red"), lty=1,lwd=1.5,bty="n")
