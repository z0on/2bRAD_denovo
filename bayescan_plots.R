source('plot_R.r')

#install.packages("boa")
library(boa)
source('plot_R.r')
dat=read.table("mydata.baye_fst.txt",header=T)
head(dat)
table(dat[,"qval"]<0.1)
outs=which(dat[,"qval"]<0.1)
plot_bayescan("mydata.baye_fst.txt",FDR=0.1,add_text=F,size=0.5,highlight=outs)

#----------------
mydata=read.table("mydata.baye.sel",colClasses="numeric")
head(mydata)
quartz()
plot(density(mydata[,parameter]), xlab=parameter, main=paste(parameter,"posterior distribution"),xlim=c(0,0.25),ylim=c(0,500))
boa.hpd(mydata[,parameter],0.05)
lines(density(mydata[,"Fst2"]))
lines(density(mydata[,"Fst3"]))
lines(density(mydata[,"Fst4"]))
boa.hpd(mydata[,parameter],0.05)

