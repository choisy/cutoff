#' @keywords internal
NULL
# This file contains utilitary functions used for the confint method of the
# S3 "em" class.

#-------------

# This function calculates the confidence intervals of parameters
# "mu1", "sigma1", "mu2" and "sigma2".
coef_ci <- function(object,level=.95) {
# "object" is an output of the function "em".
  require(bbmle) # "object$out" is of class "mle2".
  with(object,{
    ci <- cbind(coef(out),confint(out,level=level))
    ci <- as.data.frame(exp(ci))
    names(ci)[1] <- "Estimate"
    return(ci)
  })
}

#-------------

# This function computes the expectation step of the EM algorithm.
estep <- function(lambda,D1,D2,mu1,sigma1,mu2,sigma2,data,t) {
# This function basically computes the expectation step of the
# E-M algorithm (see above).
  lambda0 <- 0
  while(abs(lambda0-lambda)>t) {
    lambda0 <- lambda
    distr1 <- lambda*D1(data,mu1,sigma1)
    distr2 <- (1-lambda)*D2(data,mu2,sigma2)
    lambda <- mean(distr1/(distr1+distr2))
  }
  return(lambda)
}

#-------------

# This function returns an estimate of mean and standard deviation of lambda
# based on Oakes 1999.
lci0 <- function(object,lambda,D1,D2,data,t) {
# "object" is a named list containing the values of the mu and sigma parameters.
  with(object, {
    tmp <- estep(lambda,D1,D2,mu1,sigma1,mu2,sigma2,data,t)
    distr1 <- tmp*D1(data,mu1,sigma1)
    distr2 <- (1-tmp)*D2(data,mu2,sigma2)
    lambda <- distr1/(distr1+distr2)
    sdlambda <- sqrt(1/(sum(lambda)/(tmp^2) + sum((1-lambda))/((1-tmp)^2)))
    return(c(lambda=tmp,sd=sdlambda))
  })
}

#-------------

# This function returns an estimate of lambda and its confidence interval:
lci <- function(object,lambda,D1,D2,data,t,level=.95) {
# "object" is a named list containing the values of the mu and sigma parameters.
  out <- lci0(object,lambda,D1,D2,data,t)
  level <- (1-level)/2
  with(as.list(out),return(c(lambda,lambda+qnorm(c(level,1-level))*sd)))
}

#-------------

# This function returns an estimate of "lambda" and its confidence interval.
lambda_ci <- function(object,t=1e-64,nb=10,level=.95) {
# "object" is an output of the function "em".
# "nb" is the number of Monte Carlo simulations.
#  require(mc2d) # for "rmultinormal".
#  dHash <- c(normal=dnorm,"log-normal"=dlnorm,gamma=dgamma,weibull=dweibull)
  with(object,{
    coef <- coef(out)
    the_names <- names(coef)
    # First, draw mu1, sigma1, mu2 and sigma2 values in a multinormal distribution:
    coef <- exp(mc2d::rmultinormal(nb,coef,as.vector(vcov(out))))
    coef <- as.list(data.frame(t(coef)))
    coef <- lapply(coef,function(x){
      names(x) <- the_names
      return(as.list(x))
    })
    # Then we calculate confidence interval of lambda:
    out <- sapply(coef,function(x)lci(x,mean(lambda),
      dHash[[D1]],dHash[[D2]],data,t,level))
    out <- t(as.matrix(apply(out,1,mean)))
    level <- 100*(1-level)/2
    # Put in shape and return the output:
    colnames(out) <- c("Estimate",paste(level,"%"),paste(100-level,"%"))
    rownames(out) <- "lambda"
    return(out)
  })
}
