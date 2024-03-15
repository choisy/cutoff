#' The mode of a probability distribution.
#'
#' @keywords internal
# This function finds the mode of distribution knowing its parameters values.
findmode <- function(D,p1,p2) {
  return(switch(D,
                normal = p1,
                "log-normal" = exp(p1-p2^2),
                gamma = (p1-1)/p2,
                weibull = ifelse(p1>1,p2*((p1-1)/p1)^(1/p1),0)
  ))
}

#----------------

#' Cutoff value given parameter values.
#'
#' @importFrom stats uniroot
#' @keywords internal
# This function returns the cutoff value given parameters values.
cutoff0 <- function(mu1,sigma1,mu2,sigma2,lambda,D1,D2,distr=2,type1=.05) {
# "distr" indicates which peak we are targetting.
  interval <- c(findmode(D1,mu1,sigma1),findmode(D2,mu2,sigma2))
#  dHash <- c(normal=dnorm,"log-normal"=dlnorm,gamma=dgamma,weibull=dweibull)
  D1 <- dHash[[D1]]
  D2 <- dHash[[D2]]
  distr1 <- function(x) lambda * D1(x, mu1, sigma1)
  distr2 <- function(x) (1 - lambda) * D2(x, mu2, sigma2)
  if(distr==1) p <- function(x) distr1(x)/(distr1(x)+distr2(x))-1+type1
  else p <- function(x) distr2(x)/(distr1(x)+distr2(x))-1+type1
  return(uniroot(p,interval)$root)
}

#----------------

#' Cutoff value for a bimodal distribution.
#'
#' \code{cutoff} returns the cutoff value from a bimodal distribution, together
#'   with its confidence interval.
#'
#' From a fitted finite mixture model, we compute the probability to belong to
#' one of the two probability distributions of the finite mixture model, as a
#' function of the datum value. This probability function is used to look for the
#' cutoff value defined as the datum value for which the probability to belong to
#' a given probability distribution is equal to a type-I error. The confidence
#' interval of this cutoff value is computed by Monte Carlo simulations where in
#' each iteration the five parameter values of the finite mixture model are
#' sampled in a multinormal distribution and then the cutoff value is computed.
#' The confidence interval of the cutoff value is computed by fitting a normal
#' distribution by maximum likelihood to the Monte-Carlo-derived values of the
#' cutoff. This last step is performed by the \code{fitdistr} function of the
#' \code{MASS} package.
#'
#' @inheritParams em
#' @inheritParams confint.em
#' @param object An output from the function \code{em}.
#' @param distr Either 1 or 2, indicates which distribution belonging the Type-I
#'   error corresponds to. 1 correspond to the first distribution in \code{object}.
#' @param type1 A numerical value between 0 and 1, the value of the type-I error.
#' @return Returns a numerical vector of lenght 3, the first value being the
#'   estimated value of the cutoff and the second and thrid values being the
#'   lower and upper bound of the confidence interval of this estimate, at the
#'   level specified by parameter \code{level}.
#' @references
#'  Trang N.V., Choisy M., Nakagomi N.T., Chinh N.T.M., Doan, Y.H., Yamashiro T.,
#'    Bryant J.E., Nakagomi O. and Anh D.D. (2015) Determination of cut-off cycle
#'    threshold values in routine RT-PCR assays to assist differential diagnosis
#'    of norovirus in children hospitalized for acute gastroenteritis. Epidemiol.
#'    Infect. In press.
#' @seealso \code{\link{em}} for fitting a finite mixture model.
#' @examples
#' # Measles IgG concentration data:
#' length(measles)
#' range(measles)
#' # Plotting the data:
#' hist(measles,100,FALSE,xlab="concentration",ylab="density",ylim=c(0,.55),
#'   main=NULL,col="grey")
#' # Estimating the parameters of the finite mixture model:
#' (measles_out <- em(measles,"normal","normal"))
#' # Adding the E-M estimated finite mixture model:
#' lines(measles_out,lwd=1.5,col="red")
#' # Estimating a cutoff value from this fitted finite mixture model:
#' (cut_off <- cutoff(measles_out))
#' polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55),
#'    col=rgb(0,0,1,.2),border=NA)
#' abline(v=cut_off[-1],lty=2,col="blue")
#' abline(v=cut_off[1],col="blue")
#' @export
# This function returns cutoff value together with confidence interval from an
# output of the "em" function.
cutoff <- function(object,t=1e-64,nb=10,distr=2,type1=.05,level=.95) {
  # "object" is an output of the "em" function.
  #  require(mc2d) # for "rmultinormal".
  #  require(MASS) # for "fitdistr".
  # The dictionary:
  #  dHash <- c(normal=dnorm,"log-normal"=dlnorm,gamma=dgamma,weibull=dweibull)
  with(object,{
    coef <- out@coef
    the_names <- names(coef)
    # First we draw values for mu1, sigma1, mu2 and sigma2 in a multinormal
    # distribution:
    coef <- exp(mc2d::rmultinormal(nb,coef,as.vector(bbmle::vcov(out))))
    coef <- as.list(data.frame(t(coef)))
    coef <- lapply(coef,function(x){
      names(x) <- the_names
      return(as.list(x))}
    )
    # Then we draw random lambda value:
    out <- sapply(coef,function(x)
      lci0(x,mean(lambda),dHash[[D1]],dHash[[D2]],data,t)) # mean and sd of lambda
    lambda <- rnorm(nb,out[1,],out[2,]) # random value of lambda
    # Put all the coefficients together:
    coef <- sapply(coef,function(x)unlist(x))
    the_names <- c(rownames(coef),"lambda")
    coef <- rbind(coef,lambda) # append to other parameters
    coef <- lapply(as.data.frame(coef),function(x){
      names(x) <- the_names
      return(as.list(x))
    })
    # Call the function "cutoff0" for each combination of parameters values:
    out <- sapply(coef,function(x)
      with(x,cutoff0(mu1,sigma1,mu2,sigma2,lambda,D1,D2,distr,type1)))
    # Calculate the mean of the cutoff value, together with its confidence interval:
    out <- MASS::fitdistr(out,"normal")
    the_mean <- out$estimate["mean"]
    level <- (1-level)/2
    level <- c(level,1-level)
    ci <- the_mean+qt(level,Inf)*out$sd["mean"]
    # Put in shape and return the output:
    out <- c(the_mean,ci)
    names(out) <- c("Estimate",paste(100*level,"%"))
    return(out)
  })
}
