#Random Matrix class


eigenSpectrum<-function(input){

  call<-match.call()
  if(class(input)=="Ser"){ser<-input$series}else{
    ser<-input}
  Vars<-apply(ser,2,var)

  ser<-apply(ser,2,function(x) x/sd(x))
  Q<-nrow(ser)/ncol(ser)
  if(Q<1){stop("submit matrix as columns-variables, rows-observations, where nrow>>ncol")}
  #compute empirical covariance matrix
  E<-empCov(ser)
  #save the variances
  C<-cov2cor(E)
  #save the eigen values
  L<-eigen(C,symmetric = T)
  sigma2<-1-max(L$values)/length(L$values) #subtract the market behaviour

    #Marcenko-Pastur upper band

    require(RMTstat) #for the MP distribution

    Lmax<-sigma2*(1+1/Q+2*sqrt(1/Q))
    Lmin<-sigma2*(1+1/Q-2*sqrt(1/Q))

    x<-seq(from=Lmin,to=Lmax,length.out = length(L$values))
    theory<-data.frame(eigT=x,probT=pmp(x,var=1,svr=Q))

    type<-"noisy"
    ans<-list("call"=call,"type"="noisy","theoretical"=theory,"Vars"=Vars,"empCor"=C,"band"=c(Lmin,Lmax),"spec"=L,"Q"=Q)
  class(ans)<-"RMT"
return(ans)
  }

###### Filter RMT ######

filterRMT<-function(input){
  call<-match.call()
  if(!(class(input)=="RMT")){input<-eigenSpectrum(input)}

  n<-length(input$spec$values)
  L<-input$spec
  eig<-L$values
  Q<-input$Q

  eig[eig<=input$band[2] & eig>=input$band[1]]<-mean(eig[eig<input$band[2]])

  DD<-diag(eig,n)

  # change intro new basis - make up - new covariance matrix
  C_hat<-(L$vectors)%*%DD%*%t(L$vectors)

  #correction Brian Lee Young Rowe
  DiagM<-diag(C_hat)%o%rep(1,n)

  C_hat<-round(C_hat/(sqrt(DiagM*t(DiagM))),digits=10)

  L<-eigen(C_hat,symmetric = T,only.values=TRUE)
  sigma2<-1-max(L$values)/length(L$values) #subtract the market behaviour

  #Marcenko-Pastur upper band


  require(RMTstat) #for the MP distribution

  Lmax<-sigma2*(1+1/Q+2*sqrt(1/Q))
  Lmin<-sigma2*(1+1/Q-2*sqrt(1/Q))

  x<-seq(from=Lmin,to=Lmax,length.out = length(L$values))
  theory<-data.frame(eigT=x,probT=pmp(x,var=1,svr=Q))


  #check
  if(is.null(input$Vars)){Vars<-diag(input$cov)}

  type<-"RMT filtered"
  ans<-list("call"=call,"type"=type,"theoretical"=theory,"spec"=L,"empCor"=C_hat,"empCov"=cor2cov(C_hat,input$Vars),"Vars"=input$Vars,"band"=c(Lmin,Lmax))
  class(ans)<-"RMT"
  return(ans)
}


compMat<-function(input){
  call<-match.call(expand.dots = T)
  if(!class(input)=="Ser"){stop("Needs object of class Ser")}

  RealCov<-input$cov
# RealCor<-cov2cor(RealCov)
  #Pearson Cov
  PearCov<-empCov(input)
#  PearCor<-cov2cor(PearCov)
  #compute RMT approx of Cov and Cor
 RMTCor<-filterRMT(input)
  RMTCov<-cor2cov(RMTCor$empCor,RMTCor$Vars)
  #Wavelet covariance

  #Max distance
  mdPear<-max(abs(PearCov-RealCov))
  mdRMT<-max(abs(RMTCov-RealCov))

  #Kullback-Leibler divergence
  KLPear<-KullbackLeibler(RealCov,PearCov)
  KLRMT<-KullbackLeibler(RealCov,RMTCov)

  #Compute Frobenious norms
  cov<-list("Pearson"=frobNorm(PearCov-RealCov),"RMT"=frobNorm(RMTCov-RealCov))
  md<-list("Pearson"=mdPear,"RMT"=mdRMT)
  KL<-list("Pearson"=KLPear,"RMT"=KLRMT)
#  cor<-list("Pearson"=frobNorm(PearCor-RealCor),"RMT"=frobNorm(RMTCor$empCor-RealCor))
  ans<-list("call"=call,"Covariance"=cov,"MaxDistance"=md,"KullbackLeibler"=KL)
  class(ans)<-"compMat"
  return(ans)
}

print.compMat<-function(object){
  cat("\n Frobenious norm of Covariance matrix:\n")
  print(unlist(object$Covariance))
  cat("\n Max distance:\n")
  print(unlist(object$MaxDistance))
  cat("\n KullbackLeibler divergence:\n")
  print(unlist(object$KullbackLeibler))
}



print.RMT<-function(object){
print(paste("Random matrix called on",object$type, "matrix type",object$definite))
cat("\n Marcenko-Pastur Bound is ")
print(object$band)
cat("\n Eigenvalues lower than the bound: ")
y<-length(object$spec$values[object$spec$values<object$band[1]])

print(y)
cat("\n Eigenvalues higher than the bound: ")
yy<-length(object$spec$values[object$spec$values>object$band[2]])
print(yy)
}

plot.RMT<-function(object){
require(magrittr)

  z<-data.frame(object$spec$values,object$theoretical$eigT)

  melt(z) %>% ggplot(.,aes(value,fill=variable))+geom_density(alpha=0.5,position = "identity")+scale_fill_manual(values=c("#00C8E0", "#EE0049"))

}


empCov<-function(inp){
  if(class(inp)=="Ser"){ser<-inp$series}
    else{if(is.matrix(inp)){ser<-inp}
      else{ser<-as.matrix(inp)
      if(!is.matrix(ser)){stop("Input needs to be Ser class or something that can be converted into matrix")}}}

  if(nrow(ser)<ncol(ser)){stop("Matrix has to have more observations than variables. Maybe transpose the matrix")}
  T<-nrow(ser)
  ser<-apply(ser,2,function(x) x-mean(x))
  C<-(1/T)*t(ser)%*%ser
  return(C)
}

cor2cov<-function(mat,vars){
  sd<-diag(sqrt(as.vector(vars)),length(vars))
cov<-round(sd%*%mat%*%sd,digits=7)
return(cov)
}

frobNorm<-function(input){
  if(!is.numeric(input)){stop("Argument must be numeric")}
  if(is.matrix(input)){

    }else{if(is.vector(input)){input<-matrix(input,ncol=1)}else{stop("Input is neither matrix nor a vector")}}

  sqrt(sum(input*input))
}

KullbackLeibler<-function(mat0,mat1){
  if(class(mat0)=="Ser"){mat0<-mat0$cov$cov}
  if(class(mat0)=="myCov"){mat0<-mat0$cov}
  if(class(mat1)=="RMT"){mat1<-mat1$empCov}

  mat0<-as.matrix(mat0)
  mat1<-as.matrix(mat1)

#cat("Taking the first argument as given matrix, P, second as approximation, Q\n")
#https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Kullback.E2.80.93Leibler_divergence
  d<-0.5*(sum(diag(solve(mat1)%*%mat0))-nrow(mat0)+log(det(mat1)/det(mat0)))
 return(d)
}


compMatEmp<-function(data,by=10,sim=FALSE,plot=T){
  #one off function only

  if(sim==TRUE){
gbm<-simGBM(secs=10,Tau=250)

  emp<-sapply(seq(from=10,to=250,by = by), function(x) KullbackLeibler(gbm$cov$cov,empCov(gbm$series[1:x,])))
  rmt<-sapply(seq(from=10,to=250,by = by), function(x) KullbackLeibler(gbm$cov$cov,filterRMT(gbm$series[1:x,])))

  plot(c(seq(from=10,to=250,by=by)/10),emp,type="l", col="blue", main="EmpCov and RMT K-L dist, S&P 500",ylab="K-L distance", xlab="Q=T/N",lwd=2)
  lines(c(seq(from=10,to=250,by=by)/10),rmt,col="red",lwd=4)
return(data.frame(emp,rmt))
}else{
if(!class(data)=="Ser"){warning("Input class should be Ser")}
  reb<-seq(from=(ncol(data$series)+1),to=nrow(data$series),by=ncol(data$series))
  emp<-sapply(reb, function(x) KullbackLeibler(data$cov,empCov(data$series[1:x,])))
  rmt<-sapply(reb, function(x) KullbackLeibler(data$cov,filterRMT(data$series[1:x,])))
  Q<-round(reb/ncol(data$series))
  if(plot){
    require(magrittr)
  p<-data.frame(emp,rmt,Q) %>% melt(id.vars="Q") %>% ggplot(aes(Q,value))+geom_line(aes(colour=variable))+
     scale_y_continuous(limits=c(0,min(max(emp),max(rmt))))+ggtitle("K-L divergence")+scale_x_continuous(breaks = Q[seq(1,max(Q),5)])
  print(p)
  }
  return(data.frame(emp,rmt,Q))
}
  }

  bootComp<-function(){
  vars<-sample(1:75,sample(5:20,1))
  cov<-list("cov"=empCov(sp500.subset[,vars]))
  x<-list("series"=sp500.subset[sample(1:20,1):sample(80:200,1),vars],"cov"=cov)
  class(x)<-"Ser"
  compMat(x)$KullbackLeibler
  }





