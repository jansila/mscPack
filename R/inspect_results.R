
#statistics
#t test of returns - significance from

#difference of ES


#stable parameters

res_stablePar<-function(data){
  m<-matrix(ncol=ncol(data),nrow=4)
  n<-m
  for(i in 1:len){
    m[,i]<-stableFit(data[,i]$returns$regular)
    n[,i]<-stableFit(data[,i]$returns$denoised)
  }
  ans<-list("regular"=m,"denoised"=n)
  ans
}



#plots

res_timeES<-function(data){
  #define colours, first two red and blue, third one is orange
  cols<-c("#00C8E0","#EE0049","#FF6621")
  require(magrittr)
  len<-dim(data)[2]
  pes<-vector("list",len)

for(i in 1:len){
cat("\ni = ",i,"\n")
  a1<-as.vector(unlist(data[,i]$regular$ES))
  a2<-as.vector(unlist(data[,i]$regular$histES))
  Rebalance<-c(1:length(a1))
  dfR<-data.frame(cbind(Rebalance,a1,a2))
  names(dfR)<-c("Rebalance","ES","histES")

  a1<-as.vector(unlist(data[,i]$denoised$ES))
  a2<-as.vector(unlist(data[,i]$denoised$histES))
  dfD<-data.frame(cbind(Rebalance,a1,a2))
  names(dfD)<-c("Rebalance","ES","histES")

  a<-melt(dfR,id.vars="Rebalance")
  b<-melt(dfD,id.vars="Rebalance")
l<-max((length(Rebalance)-1),1)

  if(length(Rebalance)>1){pes[[i]]<- a %>% ggplot(.,aes(Rebalance,value))+geom_line(aes(colour=variable))+
    geom_line(data=b,aes(Rebalance,value,colour=variable),linetype=2)+
    #geom_density(data=b,aes(value,colour=variable),linetype=2)+
    scale_color_manual(values = cols[1:2])+
    scale_x_continuous(breaks=c(seq(from=1,to=Rebalance[l],by=5),Rebalance[(l+1)]))+
    ggtitle(paste("ES and realised ES, Q: ",i+1 ))
  }else{
      pes[[i]]<- a %>% ggplot(.,aes(Rebalance,value))+geom_point(aes(colour=variable))+
        geom_point(data=b,aes(Rebalance,value,colour=variable),shape=18)+
        #geom_density(data=b,aes(value,colour=variable),linetype=2)+
        scale_color_manual(values = cols[1:2])+
        ggtitle(paste("ES and realised ES, Q: ",i+1 ))
    }


}
  pes}

res_histES<-function(data){
  #define colours, first two red and blue, third one is orange
  cols<-c("#00C8E0","#EE0049","#FF6621")
  require(magrittr)
  len<-dim(data)[2]
  pes<-vector("list",len)

  for(i in 1:len){
    cat("\ni = ",i,"\n")
    a1<-as.vector(unlist(data[,i]$regular$ES))
    a2<-as.vector(unlist(data[,i]$regular$histES))
   # Rebalance<-c(1:length(a1))
    dfR<-data.frame(cbind(a1,a2))
    names(dfR)<-c("ES","histES")

    a1<-as.vector(unlist(data[,i]$denoised$ES))
    a2<-as.vector(unlist(data[,i]$denoised$histES))
    dfD<-data.frame(cbind(a1,a2))
    names(dfD)<-c("ES","histES")

    a<-melt(dfR)
    b<-melt(dfD)

  if(length(a1)>1){
    pes[[i]]<- a %>% ggplot(.,aes(value))+geom_density(aes(colour=variable))+
      geom_density(data=b,aes(value,colour=variable),linetype=2)+
      #geom_density(data=b,aes(value,colour=variable),linetype=2)+
      scale_color_manual(values = cols[1:2])+
    #  scale_x_continuous(breaks=c(seq(from=1,to=Rebalance[l],by=5),Rebalance[(l+1)]))+
      ggtitle(paste("ES and realised ES, Q: ",i+1 ))
    }else{
      pes[[i]]<- a %>% ggplot(.,aes(Rebalance,value))+geom_point(aes(colour=variable))+
        geom_point(data=b,aes(Rebalance,value,colour=variable),shape=18)+
        #geom_density(data=b,aes(value,colour=variable),linetype=2)+
        scale_color_manual(values = cols[1:2])+
        ggtitle(paste("ES and realised ES, Q: ",i+1 ))
    }


  }
  pes}

res_timeCumret<-function(data){
  #define colours, first two red and blue, third one is orange
  cols<-c("#00C8E0","#EE0049","#FF6621")
  require(magrittr)
  len<-dim(data)[2]
  pes<-vector("list",len)

  for(i in 1:len){
    l<-nrow(data[,i]$returns)
    df<-data.frame(cbind(1:nrow(data[,i]$returns),cumsum(data[,i]$returns)))
    names(df)<-c("Date","regular","denoised")
    pes[[i]]<-df %>% melt(id.vars="Date",value.name="cumret") %>%
      ggplot(.,aes(Date,cumret))+geom_line(aes(colour=variable),size=0.25)+
      scale_color_manual(values=cols[2:1])+
      ggtitle(paste("Cumulative returns, Q: ",i+1))

  }
  pes}

res_histret<-function(data){
  #define colours, first two red and blue, third one is orange
  cols<-c("#00C8E0","#EE0049","#FF6621")
  require(magrittr)
  len<-dim(data)[2]
  pes<-vector("list",len)

  for(i in 1:len){
    df<-data.frame(data[,i]$returns)
    pes[[i]]<-df %>% melt %>% ggplot(.,aes(value))+
      geom_histogram(aes(fill=variable),position="identity",alpha=0.6,bins=20)+
      scale_fill_manual(values=cols[2:1])+
      ggtitle(paste("Histogram of returns, Q: ",i+1))


  }
  pes}
