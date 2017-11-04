plotAdmixture=function(data,npops,colors=NA,ord=NULL,angle=0,space=0,grouping.method="distance",vshift=0,hshift=0,...) {
#data=tbl;npops=2;colors=NA;angle=0;ord=NULL;cluster.method="average"
	tbl=data
	require(colorspace)
	if (is.na(colors[1])){
		 colors=diverge_hcl(n=npops,h = c(260, 0), c = 100, l = c(50, 90), power = 1, fixup = TRUE, gamma = NULL, ...)
	}
	p=levels(data$pop)[1]
	tsort=c()
	if (is.null(ord)) {
		for(p in levels(tbl$pop)){
			s=subset(tbl,pop==p)
	#		dd=dist(s[,1:npops])
	#		dd[dd<0]=0
	#		dd=1-dd
			if (nrow(s)>2) {
				if (grouping.method=="distance") {
					cl=hclust(dist(s[,1:npops]))
				} else {
					if (grouping.method=="correlation") {
						dd=cor(t(s[,1:npops]))
						dd[dd<0]=0
						dd=1-dd
						cl=hclust(as.dist(dd))
					} else { stop("grouping.method not recognized\n") }
				}
				sorted=cl$labels[cl$order]
				s=s[sorted,]
			}
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
