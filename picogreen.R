# edit the  first line to match your input file name (2-column csv: well,fluorescence. No header)
# then mark and execute the whole script
# the output is printed out and saved under filename like input with prefix "res_"

inFile="picogreen.csv"

pg=read.csv(inFile,header=F) 
names(pg)=c("well","fluo")
pg[,2]=pg[,2]-pg[pg$well=="A1",2]
pg=pg[pg$well!="A1",]
pg=pg[pg$fluo>0,]
pg$fluo=log(pg$fluo,10)
cal=pg[grep("[ABCDEFGH]1$",pg$well),]
cal$conc=log(c(100,100/3,100/3^2,100/3^3,100/3^4,100/3^5,100/3^6),10)
str(cal)
plot(conc~fluo,cal,xaxt="n",yaxt="n",bty="n",xlab="fluorescence",ylab="DNA concentration (ng/ul)",mgp=c(2.1,1,0))
ll=lm(conc~fluo,cal)
abline(ll,col="red")
axis(1,at=cal$fluo,labels=round(10^(cal$fluo),0))
axis(2,at=cal$conc,labels=signif(10^(cal$conc),2))
sam=pg[-grep("[ABCDEFGH]1$",pg$well),]
sam=sam[sam$fluo>0,]
sam$conc=predict(ll,newdata=data.frame(fluo=sam$fluo))
points(conc~fluo,sam,col="cyan3")
legend("topleft",col=c("black","cyan3"),pch=1,c("calibrarion","samples"),bty="n")
sam$conc=round(10^sam$conc,1)
sam$fluo=NULL
names(sam)[2]="ng/ul"
sam
write.csv(sam,file=paste("res_",inFile,sep=""),quote=F,row.names=F)
