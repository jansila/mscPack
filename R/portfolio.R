globMin<-function(cov,e=NULL){
  call<-match.call()
  if(is.null(e)){e<-rep(0.05,dim(cov)[2])}
  #label the assets
  if(is.null(names(e))){
    asset.names <- rep("asset",length(e))
    pnumb <- c(1:length(e))
    asset.names <- paste(asset.names,pnumb)
  }else{asset.names<-names(e)}

  #inverse matrix
  cov.inv<-solve(cov)
  #solve the global min variance portfolio
  w<-as.vector(rowSums(cov.inv)/sum(cov.inv))
  names(w)<-asset.names
  eR<-crossprod(w,e)
  sd<-sqrt(t(w)%*%cov%*%w)
  globMin<-list("call"=call,"eR"=eR,"sd"=sd,"weights"=w)
  class(globMin)<-"portfolio"
  return(globMin)
}


getPortf<-function(cov,e,w){
  if(!matrixcalc::is.positive.semi.definite(cov)){stop("Covariance matrix has to be PSD")}

  n<-dim(cov)[1]
  #label the assets
  if(is.null(names(e))){
    asset.names <- rep("asset",length(e))
    pnumb <- c(1:length(e))
    asset.names <- paste(asset.names,pnumb)
  }else{asset.names<-names(e)}

    e<-matrix(e,nrow=n)
    w<-matrix(w,nrow=n)
    names(w)<-asset.names
  er<-crossprod(e,w)
  sd<-sqrt(t(w)%*%cov%*%w)
  pfolio<-list("eR"=er,"sd"=sd,"weights"=w)
  class(pfolio)<-"portfolio"
return(pfolio)
  }

effPortf<-function(cov=NULL,e=NULL,target=NULL){
  call<-match.call()
  if(is.null(cov)){cov<-simCov(20)$cov}
  if(is.null(e)){e<-(runif(nrow(cov))*0.1)}
  if(is.null(target)){target<-max(e)}
  n<-dim(cov)[1]
  e<-matrix(e,ncol=1)

  #label the assets
  if(is.null(names(e))){
    asset.names <- rep("asset",length(e))
    pnumb <- c(1:length(e))
    asset.names <- paste(asset.names,pnumb)
  }else{asset.names<-names(e)}


  unit<-rep(1,n)
  #compose matrix A then solve system
  top<-cbind(2*cov,e,unit)
  bottom<-cbind(rbind(t(e),t(unit)),matrix(0,ncol=2,nrow=2))
  A<-rbind(top,bottom)
  #prepare right handside
  b0<-matrix(c(rep(0,n),target,1),ncol=1)
  w<-solve(A,b0)[1:n]
 #compute characteristics
  er<-crossprod(e,w)
  sd<-sqrt(t(w)%*%cov%*%w)
  pfolio<-list("call"=call,"eR"=er,"sd"=sd,"weights"=w,"cov"=cov)
  class(pfolio)<-"portfolio"
  return(pfolio)
  }



tanPortf <-
  function(cov,e,risk.free)
  {

    call <- match.call()

    #
    #label the assets
    if(is.null(names(e))){
      asset.names <- rep("asset",length(e))
      pnumb <- c(1:length(e))
      asset.names <- paste(asset.names,pnumb)
    }else{asset.names<-names(e)}

    # remark: could use generalized inverse if cov.mat is positive semi-definite

    #
    # compute global minimum variance portfolio
    #
    gmin.port <- globMin(cov,e)
    if(gmin.port$eR < risk.free)
      stop("Risk-free rate cannot be greater than expected return on global minimum variance portfolio")

    #
    # compute tangency portfolio
    #
    cov.inv <- solve(cov)
    w.t <- cov.inv %*% (e - risk.free) # tangency portfolio
    w.t <- as.vector(w.t/sum(w.t))	# normalize weights
    names(w.t) <- asset.names
    er.t <- crossprod(w.t,e)
    sd.t <- sqrt(t(w.t) %*% cov %*% w.t)
    tan.port <- list("call" = call,
                     "eR" = as.vector(er.t),
                     "sd" = as.vector(sd.t),
                     "weights" = w.t)
    class(tan.port) <- "portfolio"
    tan.port
  }


effFront <-
  function(cov=NULL,e=NULL, nport=20, alpha.min=-0.1, alpha.max=1.1)
  {
  call<-match.call()
    #
  #introduce for "Ser class"
  if(is.null(cov)){cov<-simCov(20)$cov}
  if(is.null(e)){e<-(runif(nrow(cov))*0.1)}

    #
    n<-dim(cov)[1]
    e<-matrix(as.vector(e),nrow=n)

    #
    # enumarate portfolio names
    #
    port.names <- rep("port",nport)
    pnumb <- c(1:nport)
    port.names <- paste(port.names,pnumb)

    if(is.null(names(e))){
      asset.names <- rep("asset",length(e))
      pnumb <- c(1:length(e))
     asset.names <- paste(asset.names,pnumb)
    }else{asset.names<-names(e)}

    #
    # compute global minimum varian

    port.gmin <- globMin(cov,e)
    w.gmin <- port.gmin$weights

    #
    # compute efficient frontier as convex combinations of two efficient portfolios
    # 1st efficient port: global min var portfolio
    # 2nd efficient port: min var port with ER = max of ER for all assets
    #
    e.max <- max(e)
    port.max <- effPortf(cov,e,e.max)
    w.max <- port.max$weights
    #create convex combination sequence
    a <- seq(from=alpha.min,to=alpha.max,length=nport)

    #can spare transpose by using outter product
    # all(a%*%t(w.gmin)==a%o%w.gmin) is TRUE
    we.mat <- a %o% w.gmin + (1-a) %o% w.max	# rows are efficient portfolios

    er.efp <- c(we.mat %*% e)						# expected returns of efficient portfolios
    #we need var-covar matrix of efficient portofilos!
    #the diagonal are variances of portfolios.
    cov.efp <- we.mat %*% cov %*% t(we.mat)
    sd.efp <- c(sqrt(diag(cov.efp)))					# std devs of efficient frontier portfolios

    dimnames(we.mat) <- list(port.names,asset.names)

    #
    # summarize results
    #
    ans <- list("call" = call,
                "eR" = er.efp,
                "sd" = sd.efp,
                "cov"=cov,
                "weights" = we.mat)
    class(ans) <- "EF"
    ans
  }


#
# print method for portfolio object
print.portfolio <- function(x, ...)
{
  cat("Calucalated:\n")
  print(x$call, ...)
  cat("\nExpected return:    ", format(x$eR, ...), "\n")
  cat("Standard deviation: ", format(x$sd, ...), "\n")
  cat("Weights:\n")
  print(round(x$weights,4), ...)
  invisible(x)
}

#
# summary method
summary.portfolio <- function(object, risk.free=NULL, ...)
  # 	risk-free rate. If not null then
  #				compute and print Sharpe ratio
  #
{
  cat("Calculated:\n")
  print(object$call)
  cat("\nExpected return:    ", format(object$eR, ...), "\n")
  cat(  "Standard deviation: ", format(object$sd, ...), "\n")
  if(!is.null(risk.free)) {
    SharpeRatio <- (object$eR - risk.free)/object$sd
    cat("Portfolio Sharpe Ratio:       ", format(SharpeRatio), "\n")
  }
  cat("Weights:\n")
  print(round(object$weights,4), ...)
  invisible(object)
}
# hard-coded 4 digits; prefer to let user control, via ... or other argument

#
# plot weights for portfolio object


plot.portfolio<-function(object, ...){
  a<-data.frame(weight=object$weights,x=names(object$weights))
  ggplot(a,aes(x,weight))+geom_bar(stat="identity")
  }

#
# print method for Markowitz object
print.EF <- function(x, ...)
{
  cat("Calculated:\n")
  print(x$call)
  xx <- rbind(x$eR,x$sd)
  row.names(xx) <- c("ER","SD")
  colnames(xx)<-dimnames(x$weights)[[1]]
  cat("\nFrontier expected returns and standard deviations\n")
  print(round(xx,4), ...)
  invisible(x)
}

#
# summary method for Efficient frontier
summary.EF <- function(object, risk.free=NULL)
{
  call <- object$call
  asset.names <- colnames(object$weights)
  port.names <- rownames(object$weights)

  if(!is.null(risk.free)) {
    # compute efficient portfolios with a risk-free asset
    nport <- length(object$eR)
    sd.max <- object$sd[1]
    sd.e <- seq(from=0,to=sd.max,length=nport)
    names(sd.e) <- port.names

    #
    # get original er and cov.mat data from call
    er <- eval(object$call$eR)
    cov <- eval(object$call$cov)

    #
    # compute tangency portfolio
    tan.port <- tanPortf(cov,e,risk.free)
    x.t <- sd.e/tan.port$sd		# weights in tangency port
    rf <- 1 - x.t			# weights in t-bills
    er.e <- risk.free + x.t*(tan.port$eR - risk.free)
    names(er.e) <- port.names

    we.mat <- x.t %o% tan.port$weights	# rows are efficient portfolios
    dimnames(we.mat) <- list(port.names, asset.names)
    we.mat <- cbind(rf,we.mat)
  }
  else {
    er.e <- object$eR
    sd.e <- object$sd
    we.mat <- object$weights
  }
  ans <- list("call" = call,
              "eR"=er.e,
              "sd"=sd.e,
              "weights"=we.mat)
  class(ans) <- "summary.EF"
  ans
}




print.summary.EF <- function(x, ...)
{
  xx <- rbind(x$eR,x$sd)
  port.names <- names(x$eR)
  asset.names <- colnames(x$weights)
  dimnames(xx)[[1]] <- c("ER","SD")
  cat("Frontier expected returns and standard deviations\n")
  print(round(xx,4), ...)
  cat("\nPortfolio weights:\n")
  print(round(x$weights,4), ...)
  invisible(x)
}
# hard-coded 4, should let user control

#
#

plot.EF<-function(object,details=FALSE){
  cols<-c("#00C8E0", "#EE0049")

  a<-data.frame(sd=object$sd,eR=object$eR)
  if(details){
  min.p<-as.data.frame(globMin(eval(object$cov),eval(object$call$e)))
  max.p<-as.data.frame(effPortf(eval(object$call$cov),eval(object$call$e),max(eval(object$call$e))))
    p<-ggplot(a,aes(sd,eR))+geom_point()+geom_point(data=min.p,aes(sd,eR,color="Min Variance Portfolio"),color=cols[2],show.legend= TRUE,size=2.5)+geom_point(data=max.p,aes(sd,eR,color="Max eR"),color=cols[1],show.legend = T,size=2.5)
    p<-p+scale_colour_manual(name="legend", values=c("blue","red"))+ggtitle("Portfolio frontier")
  return(p)
  }
  ggplot(a,aes(sd,eR))+geom_point()
}

as.data.frame.portfolio<-function(object){
  data.frame(sd=object$sd,eR=object$eR)
}

