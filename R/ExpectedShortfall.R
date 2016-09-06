calcES<-function(input,signif=0.01){
estim<-mleFit(input)

es<-esStable(signif,alpha=estim$fit$estimate[1],beta=estim$fit$estimate[2],gamma=estim$fit$estimate[3],delta=estim$fit$estimate[4])
#if(es$abs.error<1){return(data.frame("Var"=es$VaR,"ES"=es$ES))}
#warning(paste("absolute error in integral: ",es$abs.error," on ",es$subdiv))
hist<-histES(input,signif)
#mc<-mcES(estim=estim)
#"mcES"=mc$ES,
#"mcVaR"=mc$VaR,
return(data.frame("ES"=es$ES,"VaR"=es$VaR,"histES"=hist$ES,"histVaR"=hist$VaR))
  }

esStable<-function(signif=0.01, alpha=1.82, beta=0, gamma=1, delta=0,simp=FALSE,plot=FALSE){
call<-match.call()

if(alpha<=1 || beta<(-1) || beta>1)stop("alpha has to be greater than 1, and or beta is out of unit circle")
if(abs(beta)<0.05){beta<-0}
if(abs(delta)<1e-6){delta<-0}

#  cat("\nparameters:", alpha, beta, gamma, delta,"\n")
  par<-c(alpha,beta,gamma,delta)
  names(par)<-c("alpha","beta","gamma","delta")

#special case of alpha=2 -> normal distribution

if(alpha==2){
    var<-qnorm(signif,mean=0,sd=sqrt(2))
    es1<-(1/(signif*sqrt(pi)))*exp(-0.25*var*var)

    es<-c(es1*gamma-delta)
    if(simp==TRUE){return(es)}
    ans<-list("ES"=es,"VaR"=-var,"call"=call,"par"=par)
    ans<-list("ES"=es,"int"=0, "VaR"=-var, "abs.error"=0,"subdiv"=0,"call"=call,"par"=par)

    class(ans)<-"ES"
    return(ans)
  }

var<-quantStable(signif,alpha=alpha,beta=beta,gamma=gamma,delta=delta)
B<-(sign(var)*beta)
Th<-(1/alpha)*(atan(B*tan(pi*alpha*0.5)))
if(abs(var)<1e-2){var<-0}
#special cases
if(var==0){

  if(beta==0){
#    cat("\nIm in var0,beta0\n")
    es<-(2*gamma(alpha/(alpha-1)))/pi
    if(simp==TRUE){return(es)}

    ans<-list("ES"=es, "VaR"=-var,"call"=call,"par"=par)
    class(ans)<-"ES"
    return(ans)  }

#  cat("Im in var=0 beta nonzero")

  es<-(2*gamma(alpha/(alpha-1))/(pi-2*Th))*(cos(Th)/(cos(alpha*Th)^(1/alpha)))
  if(simp==TRUE){return(es)}

  ans<-list("ES"=es, "VaR"=-var,"call"=call,"par"=par)
  class(ans)<-"ES"
  return(ans)
}

#symmetric
if(beta==0){
# cat("Im in beta 0, var nonzero case\n")

  # functions for case beta == 0
  Vb<-function(theta,alpha){((cos(theta)/sin(alpha*theta))^(alpha/(alpha-1)))*(cos(theta*(alpha-1))/cos(theta))}
  gb<-function(theta,alpha){((sin((alpha-2)*theta))/sin(alpha*theta))-(alpha*cos(theta)*cos(theta))/(sin(alpha*theta)*sin(alpha*theta))}
  ib<-function(theta,alpha,var){
    gb(theta,alpha)*exp(Vb(theta,alpha)*(-abs(var)^(alpha/(alpha-1))))
  }


  int<-integrate(ib,lower=0,upper=pi/2,alpha=alpha,var=var, rel.tol=1e-200,abs.tol=1e-200,stop.on.error = F)
#  cat("\nvalue of the integral is: ", int$val)

  es<-(alpha/(1-alpha))*(abs(var)/(pi*signif))*(int$value)*gamma-delta
  if(simp==TRUE){return(es)}

  ans<-list("ES"=es, "int"=int$value, "VaR"=-var, "abs.error"=int$abs.error,"subdiv"=int$subdivisions,"call"=call,"par"=par)
  class(ans)<-"ES"
  return(ans)

}

#general case, beta nonzero, var nonzero

integ<-function(theta,alpha,Th,var){
  ((sin(alpha*(Th+theta)-2*theta)/sin(alpha*(Th+theta)))-((alpha*cos(theta)^2)/sin(alpha*theta)^2))*
    exp((-abs(var)^(alpha/(alpha-1)))*(cos(alpha*Th)^(1/(alpha-1)))*(cos(theta)/sin(alpha*(Th+theta)))^(alpha/(alpha-1))*(cos(alpha*Th+(alpha-1)*theta)/cos(theta)) )

}

V<-function(theta,Th,alpha){
  (cos(alpha*Th)^(1/(alpha-1)))*
    ((cos(theta)/sin(alpha*(Th+theta)))^(alpha/(alpha-1)))*
    (cos(alpha*Th+(alpha-1)*theta)/cos(theta))
}

g<-function(theta,alpha,Th){
  (sin(alpha*(Th+theta)-2*theta)/sin(alpha*(Th+theta)))-
    (alpha*((cos(theta)*cos(theta)))/(sin(alpha*(Th+theta))*sin(alpha*(Th+theta))))
}

#   cat("\ngoing to integrate the general case, Th: ", round(Th,digits=4), "alpha: ",alpha, " var ",round(var,digits=4))
if(plot==TRUE){
#plot the functions
    #dev.new()
par(mfrow=c(2,2))
plot(seq(-1,1,0.1),((atan((-sign(var))*(seq(-1,1,0.1))*tan(pi*alpha*0.5)))/alpha),ylab="Theta0",type="l",lwd=3,col="green",main=paste("Theta0 of beta, alpha: ", alpha," var: ", round(var,digits=3)))
plot(seq(-Th,pi/2,length.out = 50),V(seq(-Th,pi/2,length.out = 50),Th,alpha),ylab="V",type="l",lwd=3,col="red", main=paste("Plot of V, alpha=",alpha," beta: ",beta))
plot(seq(-Th,pi/2,length.out = 50),g(seq(-Th,pi/2,length.out = 50),Th,alpha),type="l",ylab="g", lwd=3,col="black", main=paste("Plot of g, alpha=",alpha," beta: ",beta))
plot(seq(-Th,pi/2,length.out = 50),integ(seq(-Th,pi/2,length.out = 50),alpha,Th,var),type="l",ylab="integral", lwd=3,col="blue", main=paste("Plot of integrand, alpha=",alpha," beta: ",beta))
}

#cat("\nintegrating: alpha: ",alpha,", beta",beta, " signif ", signif,"\n")

ii<-integrate(integ,lower=-Th, upper=pi/2, Th=Th, alpha=alpha, var=var,
               stop.on.error = FALSE,subdivisions=200L)

es1<-((abs(var)/(signif*pi))*(alpha/(1-alpha))*ii$value)

#cat("\nintegrate value:",ii$value," error: ", ii$abs.error,"\n")

### Large integration error handle #####
if(ii$abs.error>1e-2){
  #round beta
  if(abs(beta)<0.9){B<-(sign(var)*round(beta,digits=1))}else{B<-(sign(var)*sign(beta)*0.9)}
  B<-(sign(var)*round(beta,digits=2))
  Th<-(1/alpha)*(atan(B*tan(pi*alpha*0.5)))
  ii<-integrate(integ,lower=-Th, upper=pi/2, Th=Th, alpha=alpha, var=var,
                rel.tol = 2e-308, abs.tol = 2e-308, stop.on.error = F,subdivisions=200L)
  #cat("\nintegrate value after rounding beta:",ii$value," error: ", ii$abs.error,"\n")

  #check new integration error
  if(ii$abs.error>1e-2){
    #try to go little further away from the problematic point
    d<-abs(Th-0)
    Th2<-Th-sign(Th)*0.2*d
    ii<-integrate(integ,lower=-Th, upper=pi/2, Th=Th, alpha=alpha, var=var,
                  rel.tol = 2e-308, abs.tol = 2e-308, stop.on.error = F,subdivisions=200L)
  #  cat("\nintegrate value after moving away from the left point:",ii$value," error: ", ii$abs.error,"\n")


            if(ii$abs.error>1e-2){
              #use quantile method

              ###----
              ##
              #if error is still large, call matlab

              #options(matlab="/Applications/MATLAB_R2015a.app/bin/matlab") # in case it doesnt recognise matlab command

              #write.csv(round(alpha,digits=4),file="alpha.csv",row.names = F)
              #write.csv(round(Th,digits=4),file="Th.csv",row.names = F)
              #write.csv(round(var,digits=4),file="var.csv",row.names = F)

              #  matlab.lines<-c("alpha = csvread('alpha.csv',1,0);",
              #                  "Th = csvread('Th.csv',1,0);",
              #                  "var = csvread('var.csv',1,0);",
              #                  "integ(alpha,Th,var);",
              #                  "csvwrite('ans.csv',ans);"
              #  )

              #  writeLines(matlab.lines, con="myscript.m")
              #system("/Applications/MATLAB_R2015a.app/bin/matlab -nodisplay -r \"run('myscript.m'); exit\"")
              #i<-read.csv(file="ans.csv",header=F)
              #  require(R.matlab)
              #  Matlab$startServer()
              #  matlab <- Matlab()
              #  isOpen <- open(matlab)

              #  evaluate(matlab, "ans=integ(alpha,Th,var)")
              #  setVariable(matlab, alpha=alpha,Th=Th,var=var)
              #  i<-getVariable(matlab, "ans")

              es<-qES(signif=signif,alpha=alpha,beta=beta,gamma=gamma,delta=delta)

              if(simp==TRUE){return(es)}

              ans<-list("ES"=es$ES[1],"int"=NULL, "VaR"=-var, "abs.error"=NULL,"subdiv"=NULL,"call"=call,"par"=par,"q"=TRUE)
              class(ans)<-"ES"
              return(ans)


            }
    else{
      #finish after moving away from point option
      es1<-((abs(var)/(signif*pi))*(alpha/(1-alpha))*as.vector(ii$value))
            es<-c(es1*gamma-delta)
            if(simp==TRUE){return(es)}

            ans<-list("ES"=es,"int"=ii$value, "VaR"=-var, "abs.error"=ii$abs.error,"subdiv"=NULL,"call"=call,"par"=par,"q"=NULL)
            class(ans)<-"ES"
            return(ans)}
  }else{

    #finish after rounding if error is good
    es1<-((abs(var)/(signif*pi))*(alpha/(1-alpha))*as.vector(ii$value))
    es<-c(es1*gamma-delta)
      if(simp==TRUE){return(es)}

  ans<-list("ES"=es,"int"=ii$value, "VaR"=-var, "abs.error"=ii$abs.error,"subdiv"=NULL,"call"=call,"par"=par,"q"=NULL)
  class(ans)<-"ES"
  return(ans)
      }



es<-c(es1*gamma-delta)
if(simp==TRUE){return(es)}

ans<-list("ES"=es,"int"=ii$value, "VaR"=-var, "abs.error"=ii$abs.error,"subdiv"=ii$subdivisions,"call"=call,"par"=par,"q"=NULL)
class(ans)<-"ES"
return(ans)
}

#if no large error
es1<-((abs(var)/(signif*pi))*(alpha/(1-alpha))*as.vector(ii$value))
es<-c(es1*gamma-delta)
if(simp==TRUE){return(es)}

ans<-list("ES"=es,"int"=ii$value, "VaR"=-var, "abs.error"=ii$abs.error,"subdiv"=NULL,"call"=call,"par"=par,"q"=NULL)
class(ans)<-"ES"
return(ans)
}

qES<-function(estim=NULL,l=501,signif=0.01, alpha=1.82, beta=0, gamma=1, delta=0){
  if(!is.null(estim)){
    alpha=estim$fit$estimate[1];beta=estim$fit$estimate[2]
    gamma=estim$fit$estimate[3];delta=estim$fit$estimate[4]  }
  probs<-seq(from=0,to=signif,length.out = l)[2:l]
  q<-quantStable(probs,alpha=alpha,beta=beta,gamma=gamma,delta=delta)
data.frame(ES=-mean(q))

}


histES<-function(series,signif=0.01){
  f<-floor(signif*length(series))
  c<-ceiling(signif*length(series))


  ns<-sort(series)[1:c]
  l<-length(ns)
  if(!l>1){return(data.frame(ES=-mean(ns[1]),VaR=-mean(ns[1])))}

  if(!f==c){return(data.frame(ES=-mean(mean(ns[1:l]),mean(ns[1:(l-1)])),VaR=-mean(ns[l],ns[(l-1)])))}
  data.frame(ES=-mean(ns[1:l]),VaR=-mean(ns[l]))


  }

mcES<-function(estim=NULL,l=2000,signif=0.01, alpha=1.82, beta=0, gamma=1, delta=0){
  if(!is.null(estim)){
    alpha=estim$fit$estimate[1];beta=estim$fit$estimate[2]
    gamma=estim$fit$estimate[3];delta=estim$fit$estimate[4]  }
  sim<-randStable(n=l,alpha=alpha,beta=beta,gamma=gamma,delta=delta)
histES(sim,signif=signif)
}
#other functions
print.ES<-function(obj){
  cat("\nFunction called: ")
  print(obj$call)
  cat("\nParameters:\n")
  print(obj$par)
  cat(paste("\nThe Expected shortfall calculated is:",round(obj$ES,digits=5)," with VaR: ", round(obj$VaR,digits=5)))
  cat("\nError of the integral is: ") ; print(obj$abs.error); cat(" on");print(obj$subdiv);cat(" subdivisions")
}



x<-seq(1.04,2,by=0.03)
y<-seq(-1,1,by=0.4)
#y<-c(-0.2,-0.1,-0.05,0,0.05,0.1,0.2)
meshES<-function(x,y,signif){
  ret<-sapply(y,function(y) sapply(x,function(x) esStable(signif=signif,alpha=x,beta=y,simp=T)))
  if(!class(ret)=="list"){colnames(ret)<-round(y,digits=2)
 rownames(ret)<-x}
ret}


quantMesh<-function(x,y,signif){
sapply(y,function(y) sapply(x,function(x) qES(signif=signif,alpha=x,beta=y)))
}



