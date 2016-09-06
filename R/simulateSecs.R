

#simulate a covariance matrix
simCov<-function(secs=100,neg.cor=0.3){
#create
  call<-match.call()
elm<-secs*(secs-1)/2
neg.cor<-round(neg.cor,1)
covs<-matrix(NA,ncol=secs,nrow=secs)
covs[upper.tri(covs)]<-round(sample(c(rep(-1,neg.cor*10),rep(1,(1-neg.cor)*10)),elm,replace=T)*abs(rnorm(elm)),digits=3)
covs[lower.tri(covs)]<-t(covs)[lower.tri(covs)]
diag(covs)<-abs(rnorm(secs))
c<-Matrix::nearPD(covs,corr=F,maxit = 500)$mat
ans<-list("call"=call,"cov"=matrix(c,ncol=secs))
class(ans)<-"myCov"
return(ans)
}

print.myCov<-function(object){
cat("\nCall:");print(object$call)
  n<-nrow(object$cov)
cat("\nProportion of negative relationships: ")
print(round(sum(object$cov[lower.tri(object$cov)]<0)/(0.5*(n*(n-1)))),digits=3)
}



simSeries<-function(cov=NULL,secs=100,length=500){
  call<-match.call()
if(is.null(cov)){cov<-simCov(secs)}
  if(length<=3*secs){length<-5*secs}
series<-mvrnorm(length,Sigma=cov2cor(cov$cov),mu=rep(0,dim(cov$cov)[1]))+matrix(rnorm(n=secs*length),nrow=length,ncol=secs)

ans<-list("call"=call,"cov"=cov,"series"=series)

class(ans)<-c("Ser")
return(ans)
}


plot.Ser<-function(object){
  require(dplyr)

bind_cols(object$Date,object$series) %>% melt(.,id="Date") %>% ggplot(.,aes(Date,value))+geom_line(aes(colour=variable))+theme(legend.position="none")+ggtitle("Plot of logarithmic returns")

}

myRes <- function(x) {
  ((max(x)-min(x))/(max(x)-100))*(x-min(x))+100
  }


plotSer<-function(object){

  m<-melt(object)
  ggplot(m,aes(x=m$Var1,y=m$value))+geom_line(aes(colour=as.factor(m$Var2)))+theme(legend.position="none")

  }

simGBM<-function(cov=NULL,secs=100,Tau=500,sigma=0.05,neg.cor=0.3,r=NULL){
  call<-match.call()
  if(is.null(cov)){
    cov<-simCov(secs,neg.cor)}

  secs=ncol(cov$cov)
  s<-simSeries(cov,secs=secs,length=Tau+1)

  dt<-1/(Tau+1)
  time<-seq(from=0,to=1,length.out = Tau+1)
  BM<-apply(s$series,2,function(x) cumsum(sqrt(dt)*x))
if(is.null(r)){r<-0.5*sigma*sigma}
  GBM<-apply(BM,2,function(x) 100*exp((cumsum((r-0.5*sigma*sigma)*time)+sigma*x)))

  ans<-list("call"=call,"series"=matrix(diff(log(GBM)),ncol=secs),"cov"=cov,"initPrices"=matrix(GBM[1,],ncol=secs))
  class(ans)<-c("Ser")
  return(ans)
}

makeSer<-function(data,prices=TRUE){
  call<-match.call()
  if(class(data)=="Ser"){stop("object is already class Ser")}
if(prices)
  {
    returnise<-function(x){diff(log(x))}

    ret<-apply(data[,sapply(data,is.numeric)],2,returnise)
     }else{ret<-data}

  cov<-empCov(ret)
  i<-t(data[1,sapply(data,is.numeric)])
  names(i)<-names(ret)
  ans<-list("call"=call,"series"=ret,"cov"=cov,"initPrices"=i,"Date"=data.frame(data[2:dim(data)[1],sapply(data,is.date)]))
 class(ans)<-c("Ser")
  return(ans)
}


summary.Ser<-function(object){
  cat("Series generated with:")
  print(object$call)
  cat("\n Jarque-Berra rejected in")
  print(checkN(object$series))
  cat("cases")
  cat("\nMean of returns:\n")
  print(colMeans(object$series))
  cat("\nVariance of returns:")
  print(apply(object$series,2,var))
  }

print.Ser<-function(object){
cat("Object Ser: ")
  print(str(object))
  }



illustrate<-function(secs){
  plotSer(simSeries(simCov(secs,0.2),10*secs))
}

checkN<-function(object){
  require(magrittr)
  require(fBasics)
  if(class(object)=="Ser")return(apply(object$series,2,function(x) jarqueberaTest(x)@test$p.value) %>% subset(.<0.05) %>% length()/ncol(object$series))
  apply(object,2,function(x) jarqueberaTest(x)@test$p.value) %>% subset(.<0.05) %>% length()/ncol(object)

}

is.date <- function(x) {inherits(x, 'Date')}






####################
#e<-diag(abs(10*runif(10000)))
#r<-100*rortho(10000)
#t<-round(t(r)%*%e%*%r,digits=5)   #A=Q^T%*%D%*%Q
#c<-cov2cor(t)


#test<-sapply(seq(from=100,by=500,to=10000), function(x) sum(abs(cov(mvrnorm(x,mu=rep(0,10),Sigma=t))-t)^2))

#ser<-mvrnorm(500,mu=rep(0,10),Sigma=t)
#ser
#cov(ser)

#t(chol(t))%*%chol(t)

###
#s<-Schur(cors)
#t<-pmax(s$T,0)
#cov<-round(s$Q%*%t%*%t(s$Q),digits=6)
#D<-matrix(rep(0,25),ncol=5,nrow=5)
#diag(D)<-sqrt(diag(cov))
#d<-solve(D)
#cor<-round(d%*%cov%*%d,digits=6) #can use cov2cor might be faster!

#is.symmetric.matrix(cor)
#is.positive.definite(cor)
#is.positive.semi.definite(cor)


