

zTest<-function(series,vec=NULL){
if(is.null(vec)){
    f<-mleFit(series)$fit$estimate
  alpha<-f[1];beta<-f[2]
  gamma<-f[3];delta<-f[4]
  }else{
  alpha<-vec[1];beta<-vec[2]
  gamma<-vec[3];delta<-vec[4]}


m<-mean(series)
chi<-mean(abs(series-m))
z<-m/chi*length(series)^(1-(1/alpha))
p<-1-probStable(z,alpha,beta,chi,0)
p
}




# Measure portfolios

rebalbyQ<-function(data,MaxQ=16,Maxpfolios=NULL){
  #create series of Q
    reb<-seq(from=(ncol(data$series)),to=nrow(data$series),by=ncol(data$series))[2:MaxQ]
  vecPortMat(data,reb,Maxpfolios)}

#t1<-proc.time()
#bea_realfullQ<-rebalbyQ(real,MaxQ=100)
#timeelapsed<-proc.time()-t1
#t2<-proc.time()
#bea_fulLQfullD<-rebalbyQ(fullD,MaxQ=26)
#timeeslapsed2<-proc.time()-t2

vecPortMat<-function(input,rebalance,Maxpfolios){
  bunny<-Vectorize(portMat,vectorize.args = c("rebalance"))
  bunny(input,rebalance,Maxpfolios)
}

portMat<-function(input,rebalance=50,Maxpfolios=NULL){

  if(!class(input)=="Ser"){stop("input needs to be class Ser. use function makeSer")}

train<-max(ncol(input$series)+1,rebalance)
if(!is.null(Maxpfolios)){
  #if number of pfolios could be large (low rebalance number) length.out Maxpfolios per each Q so it doens go to 100 pfolios
  rebals.end<-seq(from=(train),by=rebalance,length.out=floor((nrow(input$series)-train)/rebalance)+1)[1:Maxpfolios]}else{
    #if Q is large, there might be only little number of rebalances, so let get max out of them
    # DO NOT USE it can take ages
rebals.end<-seq(from=(train),by=rebalance,length.out=floor((nrow(input$series)-train)/rebalance)+1)
  }
#just to make sure
if(max(rebals.end)>nrow(input$series)){rebals.end<-seq(from=(train),by=rebalance,length.out=floor((nrow(input$series)-train)/rebalance)+1)}

cat("\ndobra pojdme do toho. Maxpfolios: ",Maxpfolios, "rebalance: ", rebalance,"\n")

cat("\nrebals.end", rebals.end,"\n")

reg<-lapply(rebals.end, function(x) feedPortf(input$series,x,rebalance))
den<-lapply(rebals.end, function(x) feedPortf(input$series,x,rebalance,denoise=T))
regR<-as.vector(do.call("cbind",lapply(reg,function(x) x$ret)))
denR<-as.vector(do.call("cbind",lapply(den,function(x) x$ret)))

cat("rebalance length:", length(rebals.end))
if(length(rebals.end)>1){
  #if there is more than one rebalance
rregES<-rbind(reg[[1]]$risk,t(sapply(rebals.end[2:length(rebals.end)], function(x) calcES(regR[1:x]))))
rdenES<-rbind(den[[1]]$risk,t(sapply(rebals.end[2:length(rebals.end)], function(x) calcES(denR[1:x]))))
}else{

    rregES<-reg[[1]]$risk
    rdenES<-den[[1]]$risk
}

ans<-list("regular"=rregES,"denoised"=rdenES,"returns"=data.frame(regular=regR,denoised=denR))
class(ans)<-"results"

#p<-data.frame("reg"=reg$return,"den"=den$return) %>% melt %>% ggplot(.,aes(variable,value))+geom_boxplot()
#print(p)
ans

}

#matrix
#period,series data, weights, return, cumret, ES

feedPortf<-function(data,rebals.end,rebalance,denoise=F){
 start<-rebals.end-rebalance
 per<-round((rebals.end-floor(ncol(data)))/rebalance)
   if(denoise){
    cov<-filterRMT(data[1:rebals.end,])$empCov
   } else{
    cov<-empCov(data[1:rebals.end,])
    }

  p<-globMin(cov)
  r<-p$weights%*%t(data[(start+1):rebals.end,])

 # ans<-data.frame(cbind("period"=matrix(rep(per,rebalance),ncol=1),"w"=matrix(rep(p$weights,rebalance),ncol=length(p$weights),byrow=T)
#                   ,"return"=t(r),"predSD"=matrix(rep(p$sd,rebalance),ncol=1),"SD"=matrix(rep(sd(r),rebalance),ncol=1)))

  es<-calcES(r)
  ans<-list("period"=per,"ret"=r,"lw"=max(p$weights),"weights"=p$weights,"SD"=p$sd,"risk"=es)



#colnames(ans)<-c("period",names(p$weights),"return","predSD","SD")
ans
  }





testEF<-function(input){

#sers<-input$series[,sample(1:60,15)]
sers<-input$series
rebals.end<-round(c(seq(65,to=1673,length.out = 16),1673))
#rebals.end<-c(round(seq(16,to=1673,by = 16)),1673)
#rebals.end<-c(round(seq(from=11,to=200, length.out = 16)),200)
ef<-vector("list",16)
efRMT<-vector("list",16)
efR<-vector("list",16)
efRMTR<-vector("list",16)
x<-rebals.end

for(i in 1:16){
ef[[i]]<-effFront(empCov(sers[1:x[i],]),colMeans(sers[1:x[i],]))
efRMT[[i]]<-effFront(eval(filterRMT(sers[1:x[i],])$empCov),as.vector(colMeans(sers[1:x[i],])))

efR[[i]]<-effFront(empCov(sers[1:x[i+1],]),colMeans(sers[1:x[i+1],]))
efRMTR[[i]]<-effFront(eval(filterRMT(sers[1:x[i+1],])$empCov),as.vector(colMeans(sers[1:x[i+1],])))
}

ans<-list("reg"=ef,"rmt"=efRMT,"regR"=efR,"rmtR"=efRMTR)
ans
}


Lal<-function(input){
  sers<-input$series
  #rebals.end<-round(seq(14,to=200,length.out = 16))
    rebals.end<-c(102,205)
  ef<-vector("list",1)
  efRMT<-vector("list",1)
  efR<-vector("list",1)
  efRMTR<-vector("list",1)
  x<-rebals.end
  for(i in 1:1){
    ef[[i]]<-effFront(empCov(sers[1:x[i],]),colMeans(sers[1:x[i+1],]))
    efRMT[[i]]<-effFront(eval(filterRMT(sers[1:x[i],])$empCov),as.vector(colMeans(sers[1:x[i+1],])))

    efR[[i]]<-effFront(empCov(sers[1:x[i+1],]),colMeans(sers[1:x[i+1],]))
    efRMTR[[i]]<-effFront(eval(filterRMT(sers[1:x[i+1],])$empCov),as.vector(colMeans(sers[1:x[i+1],])))
    cat("\n i:",i)
  }

  ans<-list("reg"=ef,"rmt"=efRMT,"regR"=efR,"rmtR"=efRMTR)
  ans
}



pl<-function(x){
plot(ef[[x]]$sd,ef[[x]]$eR, xlim=c(min(ef[[x]]$sd,efRMT[[x]]$sd),max(ef[[x]]$sd,efRMT[[x]]$sd)),
     ylim=c(min(ef[[x]]$eR,efRMT[[x]]$eR),max(ef[[x]]$eR,efRMT[[x]]$eR)),col="blue",main=paste())
points(efRMT[[x]]$sd,efRMT[[x]]$eR,col="red")}


plotGG<-function(){

  plots<-vector("list",length(rebals.end))
  for(x in 1:length(rebals.end)){

efF<-data.frame(sd=ef[[x]]$sd,eR=ef[[x]]$eR)
efFRMT<-data.frame(sd=efRMT[[x]]$sd,eR=efRMT[[x]]$eR)
p<-ggplot(efF,aes(sd,eR,frame=rebal))+geom_point(aes(colour="EmpCov"))+geom_point(data=efFRMT,aes(sd,eR,colour="RMT"))+scale_fill_manual(values=c("#00C8E0", "#EE0049"))
plots[[x]]<-p
rm(p)
  }
return(plots)
}






myacf<-function(series,k=30){
  call<-match.call()
  bea_is_seksi<-as.vector(acf(series, lag.max=k, plot=FALSE)$acf)
  x<-seq(0,length(bea_is_seksi)-1,1)
  bunz<-cbind(x,bea_is_seksi)
  lim1<-min(-1.96/sqrt(length(series)),min(bea_is_seksi))
  lim2<-max(1.96/sqrt(length(series)),max(bea_is_seksi))
  plot(bunz, ylim=c(lim1,lim2), type="h",col="darkblue",lwd=3, xlab="Lag",ylab="ACF", main=call)
  abline(h=0, lwd=2)
  abline(h=1.96/sqrt(length(series)), col="darkred", lwd=2)
  abline(h=-1.96/sqrt(length(series)), col="darkred", lwd=2)
  points(bunz, pch=16, col="darkblue",lwd=3)
}

mypacf<-function(series,k=30){
  bea_is_seksi<-as.vector(pacf(series, lag.max=k, plot=FALSE)$acf)
  x<-seq(0,length(bea_is_seksi)-1,1)
  bunz<-cbind(x,bea_is_seksi)
  lim1<-min(-1.96/sqrt(length(series)),min(bea_is_seksi))
  lim2<-max(1.96/sqrt(length(series)),max(bea_is_seksi))
  plot(bunz, type="h",ylim=c(lim1,lim2),col="darkblue",lwd=3, xlab="Lag",ylab="PACF")
  abline(h=0, lwd=2)
  abline(h=1.96/sqrt(length(series)), col="darkred", lwd=2)
  abline(h=-1.96/sqrt(length(series)), col="darkred", lwd=2)
  points(bunz, pch=16, col="darkblue",lwd=3)
}




# save in file

#dev.off()
#pdf("plot_EFF_points")
#par(mfrow=c(2,2))
#sapply(c(1:length(rebals.end)), pl)
#dev.off()
#pdf(file="ggplot_EFF", paper="a4r")
#myPack::multiplot(plotlist = plotGG(), cols = 5)


# sort out the variance structure to get RMT covariance!



