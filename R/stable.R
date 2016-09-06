dyn.load("src/stable.so")


###########################################################################
#' Density of a stable distribution
#'
#' This function produces the stable density computed at y.
#' Note that it is obtained by integrating a function from 0 to Infinity.
#' This integral is approximated by the finite integral from 0 to UP
#' using the Simpson's method with npt points or Romberg's integration. Core of the code is based on Philipe Lambert's code https://github.com/cran/stable
#'
#'
#' @param input series of input
#' @param vec vector of class "stableFit" MLE estimates from mleFit
#' @return vector of densities of input for given 4 parameters (alpha-tail, beta-skew, gamma-scale, delta-location)
#' @examples denStable(rnorm(500), alpha=2,beta=0,gamma=1/sqrt(2), delta=0)
#'


#' @useDynLib stable stable
denStable <- function(input, vec=NULL, alpha=1.8, beta=0, gamma=1, delta=0,
                    npt=501, up=10, eps=1.0e-6){
  if(!is.null(vec) && class(vec)=="stableFit" && !any(is.na(vec$estimate)) ){
    alpha=c(vec$estimate[1])
    beta=c(vec$estimate[2])
    gamma=c(vec$estimate[3])
    delta=c(vec$estimate[4])
    }else{alpha=alpha;beta=beta;gamma=gamma;delta=delta}


  if(any(is.na(c(alpha,beta,gamma,delta)))){
        alpha=1.8;beta=0;gamma=1;delta=0}

#special cases
#  if (alpha == 2 && beta==0) {
#        return(stats::dnorm(y, mean = 0, sd = 1))
#  } else{ if (alpha == 1 && beta == 0) {
#        return(stats::dcauchy(y, log=log))}}


  ly <- length(input)
  z0 <- rep(0,ly)
  skew <- beta+z0
  tail <- alpha+z0
  yy <- (input-delta)/gamma
  z <- .C("stable",
	  as.integer(ly),
	  yy,
	  skew,
	  tail,
	  as.integer(npt),
	  up,
	  eps,
	  err=integer(1),
	  ffy=double(ly))

  z$ffy/gamma
  }


###########################################################################

#' Function computes


#' @useDynLib stable stable
probStable <- function(y, loc=0,disp=1,skew=0,tail=2,eps=1.0e-6){
  if (alpha == 2 && beta==0) {
    return(pnorm(y, mean = 0, sd = 1, log=log))
  } else if (alpha == 1 && beta == 0) {
    return(pcauchy(y, log=log))}



  yy <- (y-loc)/disp
  ly <- length(yy)
  z0 <- rep(0,ly)
  skew <- skew+z0
  tail <- tail+z0
  eta <- skew*(1-abs(1-tail))*pi/2
       if ((yy==0)&(eta==0)) return(0.5+z0)

    z <- .C("pstable",
	  as.integer(ly),
	  yy,
	  skew,
	  tail,
	  eps,
	  err=integer(1),
	  ffy=double(ly))
  z$ffy}

###########################################################################
# Quantile function of a stable random variable
#

quantStable <- function(p, alpha=1.8,beta=0,gamma=1,delta=0){
 if(alpha==2){return(qnorm(p))}

  if(alpha==1 && beta){return(qcauchy(p))}

  q<-Vectorize(stabledist::qstable,vectorize.args = c("p","alpha","beta","gamma","delta"))

    q(p,alpha,beta,gamma,delta,pm=1)

  }

###########################################################################
# Generation of stable random deviates
#

randStable <- function(n=1,alpha=2,beta=0,delta=0,gamma=1){
 return(quantStable(runif(n),delta=delta,gamma=gamma,beta=beta,alpha=alpha))}
