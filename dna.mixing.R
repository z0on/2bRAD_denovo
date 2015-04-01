dna.mixing=function(model,data,prob=0.75,...) {
	require(MCMCglmm)
	allmean=apply(model$Sol,2,mean)
	sams=grep("sam.\\(Intercept\\).sam.",names(allmean))
	pattern="sam.\\(Intercept\\).sam."
	if (is.na(sams[1])) {
		sams=grep("sam.",names(allmean))
		pattern="sam."
	}
	means=allmean[sams]
	nsam=length(means)
	hpds=HPDinterval(model$Sol[,sams],prob=prob)
	row.names(hpds)=sub(pattern,"",row.names(hpds))
	names(means)=sub(pattern,"",names(means))
	gmean=ggr=gsam=c()
	ghpd=c()
	for (g in levels(data$lane)) {
		sub=subset(data,lane==g)
		gn=levels(as.factor(as.character(as.vector(sub$sam))))
		nn=length(gn)
		meann=means[gn]
		mnn=min(meann)
		hpdd=hpds[gn,]
		meann=(meann-mnn)
		hpdd[,1]=(hpdd[,1]-mnn)
		hpdd[,2]=(hpdd[,2]-mnn)
		gmean=append(gmean,meann)
		ghpd=rbind(ghpd,hpdd)
		ggr=append(ggr,rep(g,nn))
		gsam=append(gsam,names(meann))
		
		gn %in% row.names(hpds)
		
		ymin=min(c(hpdd[,1],hpdd[,2]))-0.5
		ymax=round(max(c(hpdd[,1],hpdd[,2])),0)+0.5
		ylabs=2^(seq(0,ymax,1))
		plot(meann~c(1:nn),xaxt="n",yaxt="n",mgp=c(2.3,1,0),ylim=c(ymin,ymax),xlab="",main=paste("lane",g),ylab="microliters to mix")
		abline(h=seq(-4,14,1),lty=3,col="grey70")
		arrows(c(1:nn),meann,c(1:nn),hpdd[,1],angle=90,code=2,length=0.03)
		arrows(c(1:nn),meann,c(1:nn),hpdd[,2],angle=90,code=2,length=0.03)
		axis(side=1,labels=names(meann),at=c(1:nn),las=2)
		axis(side=2,labels=ylabs,at=seq(0,ymax,1),las=2)
	}
	res=data.frame(cbind("ul2mix"=2^gmean,2^ghpd))
	res$lane=ggr
	return(res)
}

#----------------------------
