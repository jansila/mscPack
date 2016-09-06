



dick2<-portMat(real,rebalance=100)

#replicate my own sp500
#sp<-makeSer(sp500.subset,prices = F)
t<-testEF(fullD)

p<-vector("list",16)
#lims<-findlim(t)
cols<-c("#00C8E0", "#EE0049")

for(i in 1:16){
  a<-data.frame(sd=t$reg[[i]]$sd,eR=t$reg[[i]]$eR)
  b<-data.frame(sd=t$rmt[[i]]$sd,eR=t$rmt[[i]]$eR)
  c<-data.frame(sd=t$regR[[i]]$sd,eR=t$regR[[i]]$eR)
  d<-data.frame(sd=t$rmtR[[i]]$sd,eR=t$rmtR[[i]]$eR)


  p[[i]]<-ggplot(a,aes(sd,eR))+geom_point(colour=cols[2])+geom_point(data = b,aes(sd,eR),colour=cols[1])+
     #lims(x=lims$xlim,y=lims$ylim)+
    geom_point(data=c,aes(sd,eR),shape=18,colour=cols[2])+geom_point(data=d,aes(sd,eR),shape=18,colour=cols[1])+
    ggtitle(paste("Q=",i))

}
multiplot(plotlist=p,cols=4)



#replicate Laloux
tt<-Lal(sim)
pp<-vector("list",1)
for(i in 1:1){
  a<-data.frame(sd=tt$reg[[i]]$sd,eR=tt$reg[[i]]$eR)
  b<-data.frame(sd=tt$rmt[[i]]$sd,eR=tt$rmt[[i]]$eR)
  c<-data.frame(sd=tt$regR[[i]]$sd,eR=tt$regR[[i]]$eR)
  d<-data.frame(sd=tt$rmtR[[i]]$sd,eR=tt$rmtR[[i]]$eR)
  pp[[i]]<-ggplot(a,aes(sd,eR))+geom_point(colour=cols[2])+geom_point(data = b,aes(sd,eR),colour=cols[1])+
    #  lims(x=lims$xlim,y=lims$ylim)+
    geom_point(data=c,aes(sd,eR),shape=18,colour=cols[2])+geom_point(data=d,aes(sd,eR),shape=18,colour=cols[1])+
    ggtitle("Replicating Laloux results")
}
pp

#
t<-testEF(sim)
p<-vector("list",16)
#lims<-findlim(t)
cols<-c("#00C8E0", "#EE0049")

for(i in 1:16){
  a<-data.frame(sd=t$reg[[i]]$sd,eR=t$reg[[i]]$eR)
  b<-data.frame(sd=t$rmt[[i]]$sd,eR=t$rmt[[i]]$eR)
  c<-data.frame(sd=t$regR[[i]]$sd,eR=t$regR[[i]]$eR)
  d<-data.frame(sd=t$rmtR[[i]]$sd,eR=t$rmtR[[i]]$eR)


  p[[i]]<-ggplot(a,aes(sd,eR))+geom_point(colour=cols[2])+geom_point(data = b,aes(sd,eR),colour=cols[1])+
 #  lims(x=lims$xlim,y=lims$ylim)+
    geom_point(data=c,aes(sd,eR),shape=18,colour=cols[2])+geom_point(data=d,aes(sd,eR),shape=18,colour=cols[1])+
    ggtitle(paste("Q=",i))

}
multiplot(plotlist=p,cols=4)


findlim<-function(t){
xlim<-c(min(sapply(t$reg, function(x) min(x$sd)),sapply(t$rmt, function(x) min(x$sd))),
        max(sapply(t$reg, function(x) max(x$sd)),sapply(t$rmt, function(x) max(x$sd))))

ylim<-c(min(sapply(t$reg, function(x) min(x$eR)),sapply(t$rmt, function(x) min(x$eR))),
        max(sapply(t$reg, function(x) max(x$eR)),sapply(t$rmt, function(x) max(x$eR))))

xxlim<-c(min(sapply(t$regR, function(x) min(x$sd)),sapply(t$rmtR, function(x) min(x$sd))),
        max(sapply(t$regR, function(x) max(x$sd)),sapply(t$rmtR, function(x) max(x$sd))))

yylim<-c(min(sapply(t$regR, function(x) min(x$eR)),sapply(t$rmtR, function(x) min(x$eR))),
        max(sapply(t$regR, function(x) max(x$eR)),sapply(t$rmtR, function(x) max(x$eR))))

xlim<-c(min(rbind(xlim,xxlim)[,1]),max(rbind(xlim,xxlim)[,2]))

ylim<-c(min(rbind(ylim,yylim)[,1]),max(rbind(ylim,yylim)[,2]))

data.frame(cbind(xlim,ylim))
}


pdf("price_check2.pdf",paper="a4r")
for(i in 2:65){
  print(ggplot(real,aes_string("Date",names(real)[i]))+geom_line())
}

dev.off()






ef<-lapply(rebals.end, function(x) effFront(empCov(sers[1:x,]),colMeans(sers[1:x,])))

efRMT<-lapply(rebals.end, function(x) effFront(eval(filterRMT(sers[1:x,])$empCov),as.vector(colMeans(sers[1:x,])))  )


glM<-lapply(rebals.end, function(x) globMin(empCov(sers[1:x,])))
glR<-lapply(rebals.end, function(x) globMin(eval(filterRMT(sers[1:x,])$empCov)))

emp<-sapply(1:length(glM),function(x) glM[[x]]$sd)
rmt<-sapply(1:length(glM),function(x) glR[[x]]$sd)
d<-data.frame(emp=emp,rmt=rmt)

melt(d) %>% ggplot(.,aes(value))+geom_density(aes(fill=variable,position="identity"),alpha=0.6)

plot(sapply(1:length(glM),function(x) glM[[x]]$sd),ylim=c(0,0.03))
dpoints(sapply(1:4,function(x) glR[[x]]$sd),col="red")


w<-matrix(glR[[1]]$weights,nrow=1)

ret<-sort(matrix(sers,ncol=10)%*%t(w))
qplot(ret)
quantile(ret,probs=0.05)
mean(ret[ret<quantile(ret,0.05)])
}
