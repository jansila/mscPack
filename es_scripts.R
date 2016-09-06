

# AMD seems quite large - checking the stable fit from fStable
#  Model:
#  Stable Distribution
#
#Estimated Parameter(s):
#  alpha          beta         gamma         delta
#1.636114e+00  7.577797e-03  1.787863e-02 -5.926047e-05



# ES and VaR for all individual stocks
desc<-apply(real$series,2,calcES)

#mle stable fits
descMLE<-apply(real$series,2,mleFit)
d<-data.frame(do.call("rbind",lapply(1:length(descMLE), function(x) descMLE[[x]]$fit$estimate)))
names(d)<-c("alpha","beta","gamma","delta")
row.names(d)<-ticks
write(stargazer::stargazer(d,summary=F),file=paste(tables,"desc_mle.tex",sep=""))
tables=c("~/Desktop/Dissertation/athesis/tables/")
