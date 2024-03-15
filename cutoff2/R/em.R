#' Minus log-likelihood of the finite mixture model.
#'
#' @keywords internal
# This function returns minus the log-likelihood of the set of parameters "mu1",
# "sigma1", "mu2" and "sigma2", given dataset "data", distributions "D1" and
# "D2", and parameter "lambda".
mLL <- function(mu1,sigma1,mu2,sigma2,lambda,data,D1,D2,penaltyScale) {
# "mu1", "mu2", "sigma1" and "sigma2" are the log of the parameter values.
# "lambda" and "data" are two numeric vectors of the same length.
# "D1" and "D2" are two functions of probability density.
  params <- c(mu1,sigma1,mu2,sigma2)
  names(params) <- c("mu1","sigma1","mu2","sigma2")
  with(as.list(exp(params)), {
    out <- lambda*D1(data,mu1,sigma1)
    out <- out+(1-lambda)*D2(data,mu2,sigma2)
    ## Penalties
    ## if (penaltyScale > 0) {
    ##   out <- out + penalty()
    ## }
    return(-sum(log(out)))
  })
}

##penalty

#-------------

#' Starting values of the parameter of a finite mixture model.
#'
#' @keywords internal
# This function calculates the starting values of parameter "lambda".
startval <- function(data,D1,D2) {
  #  require(tree) # for "tree".
  thresh <- tree::tree(data~data)$frame$yval[1]
  sel <- data<thresh
  data1 <- data[sel]
  data2 <- data[!sel]
  lambda <- length(data1)/length(data)
  param1 <- MASS::fitdistr(data1,D1)$est
  param2 <- MASS::fitdistr(data2,D2)$est
  out <- c(param1,param2,lambda)
  names(out) <- c("mu1","sigma1","mu2","sigma2","lambda")
  return(out)
}

#-------------

#' Expectation-Maximization estimation of a finite mixture model
#'
#' \code{em} returns points estimations of the parameters of a finite mixture
#'		model using the Expectation-Maximization (E-M) algorithm.
#'
#' The finite mixture model considered in this function is a mixture of two
#'  probability distributions that are one of the following: normal, log-normal,
#'  gamma or Weibull. Each of these distributions is defined by two parameters:
#'  a location and a scale parameter:
#' \tabular{lcc}{
#'              \tab location \tab scale \cr
#'   normal     \tab mean     \tab sd    \cr
#'   log-normal \tab meanlog  \tab sdlog \cr
#'   gamma      \tab shape    \tab rate  \cr
#'   Weibull    \tab shape    \tab scale
#' }
#' These parameters, together with the mixture parameter, are estimated by the
#'  Expection-Maximization algorithm.
#'
#' @param data A vector of real numbers, the data to model with a finite mixture
#'  model.
#' @param D1,D2 Probability distributions constituting the finite mixture model.
#'  See Details.
#' @param t A numerical scalar indicating the value below which the E-M
#' 	algorithm should stop.
#' @return A list with class \code{em} containing the following components:
#' 	\item{lambda}{a numerical vector of length \code{length(data)} containing,
#'    for each datum, the probability to belong to distribution \code{D1}.}
#'  \item{param}{the location (mu) and scale (sigma) parameters of the
#'	  probability distributions \code{D1} and \code{D2}.}
#'    \code{mu2} and \code{lambda2}.}
#' 	\item{D1,D2}{character scalars containing the names of the probability
#'    distributions used in the finite mixture model.}
#' 	\item{data}{the numerical vector of data used as input.}
#' 	\item{data_name}{character scalar containing the name of the dataset used as
#'    input.}
#' 	\item{out}{an object of class \code{mle2} that contains the
#'    maximum-likelihood estimates of parameters \code{mu1}, \code{lambda1},
#' 	\item{t}{the input \code{t} argument value.}
#' @references
#' 	Chuong B. Do and Serafim Batzoglou (2008) What is the expectation
#'    maximization algorithm? Nature Biotechnology 26(8): 897-899.\cr
#'  \cr
#'  Peter Schlattmann (2009) Medical Applications of Finite Mixture Models.
#'    Springer-Verlag, Berlin.
#' @seealso \code{\link{confint.em}} method for calculating the confidence
#'    intervals of the parameters and \code{\link{cutoff}} for deriving a
#'    cut-off value.
#' @examples
#' # Measles IgG concentration data:
#' length(measles)
#' range(measles)
#' # Plotting the data:
#' hist(measles,100,F,xlab="concentration",ylab="density",ylim=c(0,.55),
#'   main=NULL,col="grey")
#' # The kernel density:
#' lines(density(measles),lwd=1.5,col="blue")
#' # Estimating the parameters of the finite mixture model:
#' (measles_out <- em(measles,"normal","normal"))
#' # The confidence interval of the parameter estimates:
#' confint(measles_out,t=1e-64,nb=100,level=.95)
#' # Adding the E-M estimated finite mixture model:
#' lines(measles_out,lwd=1.5,col="red")
#' # The legend:
#' legend("topleft",leg=c("non-parametric","E-M"),col=c("blue","red"),
#'   lty=1,lwd=1.5,bty="n")
#' @export
# This function uses the EM algorithm to calculates parameters "lambda"
# (E step), "mu1", "sigma1", "mu2" and "sigma2" (M step).
em <- function(data, D1, D2, t=1e-64, penaltyScale=0) {
  data_name <- unlist(strsplit(deparse(match.call()),"="))[2]
  data_name <- sub(",.*$","",gsub(" ","",data_name))
  start <- as.list(startval(data,D1,D2))
  D1b <- dHash[[D1]]
  D2b <- dHash[[D2]]
  lambda0 <- 0 # the previous value of lambda (scalar).
  with(start, {
    while(abs(lambda0-mean(lambda))>t) {
      lambda <- mean(lambda)
      lambda0 <- lambda
# Expectation step:
      distr1 <- lambda*D1b(data,mu1,sigma1)
      distr2 <- (1-lambda)*D2b(data,mu2,sigma2)
      lambda <- distr1/(distr1+distr2) # lambda is a vector.
# Minimization step (maximum-likelihood parameters estimations):
      mLL2 <- function(mu1,sigma1,mu2,sigma2)
			return(mLL(mu1,sigma1,mu2,sigma2,lambda,data,D1b,D2b,penaltyScale))
      start <- as.list(log(c(mu1=mu1,sigma1=sigma1,mu2=mu2,sigma2=sigma2)))
      out <- bbmle::mle2(mLL2,start,"Nelder-Mead")
# The following 4 lines assign the MLE values to the corresponding parameters:
      coef <- out@coef
      coef_n <- names(coef)
      names(coef) <- NULL
      for(i in 1:4) assign(coef_n[i],exp(coef[i]))
    }
# Put in shape and return the output:
    out <- list(
		lambda=lambda,param=exp(out@coef),D1=D1,D2=D2,deviance=out@min,
			data=data,data_name=data_name,out=out,t=t)
    class(out) <- "em"
    return(out)
  })
}

#-------------

#' Print method of S3-class "em".
#'
#' @export
#' @method print em
print.em <- function(object) {
  hash <- list(
    normal=c("mean","sd"),
    "log-normal"=c("meanlog","sdlog"),
    gamma=c("shape","rate"),
    weibull=c("shape","scale")
  )
  param <- as.list(object$param)
  digits <- unlist(options("digits"))
  cat(paste0("\n        Finite mixture model fitting to dataset \"",
		object$data_name,"\":\n\n"))
  cat(paste("probability to belong to distribution 1 =",
		round(mean(object$lambda),digits),"\n"))
  cat(paste("distribution 1:",object$D1,"\n"))
  cat(paste0("  location (",hash[[object$D1]][1],") = ",
		round(param$mu1,digits),"\n"))
  cat(paste0("  scale (",hash[[object$D1]][2],") = ",
		round(param$sigma1,digits),"\n"))
  cat(paste("distribution 2:",object$D2,"\n"))
  cat(paste0("  location (",hash[[object$D2]][1],") = ",
		round(param$mu2,digits),"\n"))
  cat(paste0("  scale (",hash[[object$D2]][2],") = ",
		round(param$sigma2,digits),"\n"))
  cat(paste("deviance of the fitted model:",
        round(object$deviance,digits),"\n\n"))
}

#-------------

#' Lines method of S3-class "em".
#'
#' @export
#' @method lines em
lines.em <- function(object,...) {
# ...: parameter passed to the "line" function.
  with(object,with(as.list(param), {
    lambda <- mean(lambda)
    curve(dnorm(x,mu1,sigma1),add=T,n=512,lty=2,...)
    curve(dnorm(x,mu2,sigma2),add=T,n=512,lty=2,...)
    curve(lambda*dnorm(x,mu1,sigma1)+
      (1-lambda)*dnorm(x,mu2,sigma2),add=T,n=512,...)
  }))
}

#-------------

#' Confint method of S3-class \code{em}.
#'
#' This method of the S3 \code{em} class calculates the confidence intervals of
#' the parameters of the fitted finite mixture model.
#'
#' Confidence intervals of the parameters of probability distributions are
#' calculated by the \code{confint} method of the S4 \code{confint} class of the
#' \code{bbmle} package with default values. See the help of this method for
#' technical details. The confidence interval of the \code{lambda} parameter is
#' calculated thanks to the information-based method of Oakes (1999). The
#' possible non-independance between lambda and the parameters of the probability
#' distributions is accounted for by Monte Carlo simulations where each iteration
#' consists in (i) sampling values of theses parameters in a multinormal
#' distribution, and (ii) applying the method of Oakes (1999). Samplings in the
#' multinormal distribution is performed by the \code{rmultinormal} function of
#' the \code{mc2d} package.
#' @inheritParams em
#' @param nb Number of Monte Carlo simulations.
#' @param level The confidence level required.
#' @return A dataframe containing point estimates (first column) and confidence
#'   intervals (second and third columns) of each of the five parameters of the
#'   finite mixture model (five row, one per parameter).
#' @references David Oakes (1999) Direct calculation of the information matrix
#'   via the EM algorithm. J R Statist Soc B, 61: 479-482.
#' @export
#' @method confint em
# This function returns the parameter values and their confidence
# intervals from an output of the "em" function.
confint.em <- function(object,t=1e-64,nb=10,level=.95) {
  a <- coef_ci(object,level)
  b <- lambda_ci(object,t,nb,level)
  return(rbind(a,b))
}
