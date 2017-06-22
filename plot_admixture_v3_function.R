plotAdmixture=function(data,npops,colors=NA,ord=NULL,angle=0,space=0,cluster.method="single",vshift=0,hshift=0,...) {
#data=tbl;npops=4;colors=cols;angle=0;ord=NULL
	tbl=data
	require(RColorBrewer)
	if (is.na(colors[1])){
		if (npops<3) {  colors=brewer.pal(n=npops+1,name="Set3") } else { colors=brewer.pal(n=npops,name="Set3") }
	}
	p=levels(data$pop)[1]
	tsort=c()
	if (is.null(ord)) {
		for(p in levels(tbl$pop)){
			s=subset(tbl,pop==p)
			dd=cor(t(s[,1:npops]))
			dd[dd<0]=0
			dd=1-dd
			cl=hclust(as.dist(dd),method=cluster.method)
			sorted=cl$labels[cl$order]
			s=s[sorted,]
			tsort=data.frame(rbind(tsort,s))
		}
	} else {
		tsort=tbl[ord,]
	}
	midpts=barplot(t(as.matrix(tsort[1:npops])), col=colors,xlab="population", space=space,ylab="Ancestry", border=NA, xaxt="n",mgp=c(2,1,0))
	pops=levels(tbl$pop)
	np=0;lp=0;lpoints=c()
	abline(v=0)
	abline(h=1,lwd=1)
	abline(h=0,lwd=1)
	for( p in pops) {
		np0=table(tbl$pop==p)[2]
		lp0=np0/2
		lp=np+lp0
		np=np+np0
		lpoints=append(lpoints,lp)
		abline(v=np)
	}
	text(as.numeric(lpoints)+3.5+hshift,par("usr")[3]-0.1-vshift,labels=pops,xpd=TRUE, srt=angle, pos=2)
	row.names(tsort)=tsort$ind
#	return(tsort)
	return(row.names(tsort))
}
