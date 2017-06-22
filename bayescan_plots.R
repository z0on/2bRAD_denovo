setwd("~/Documents/RAD_2digitifera_4dadi/apr27_unifiedGenotyper/genomeScans/")

source('plot_R.r')

#install.packages("boa")
library(boa)


dat=read.table("snp.baye_fst.txt",header=T)
head(dat)
table(dat[,"qval"]<0.1)
# FALSE  TRUE 
# 24736   119 
# map-anchored: 9244    32 
outs=which(dat[,"qval"]<0.1)
plot_bayescan("snp.baye_fst.txt",FDR=0.1,add_text=F,size=0.5,highlight=outs)

toRemove=which(dat$qval<0.1)
toRemove

# splitting the vcf file into header and body 
system("cat gatk2bs.recode.vcf | grep '^#' >header.vhap")
system("cat gatk2bs.recode.vcf | grep -v '^#' >body.vhap")

# reading body, removing rows
vcf=read.table("body.vhap",header=F,sep="\t")

nrow(vcf)
vcf=vcf[-toRemove,]
nrow(vcf)

# writing the new body (with rows removed) down
write.table(vcf,file="body1.vhap",quote=F,row.names=F,sep="\t")
# removing header line
system("tail -n +2 body1.vhap >nohead.vhap")
# combining header and new body
system("cat header.vhap nohead.vhap >noOuts.hap.vcf")

#----------------
mydata=read.table("snp.baye.sel",colClasses="numeric")
head(mydata)
parameter="Fst1"
quartz()
plot(density(mydata[,parameter]), xlab=parameter, main=paste(parameter,"posterior distribution"))
boa.hpd(mydata[,parameter],0.05)

