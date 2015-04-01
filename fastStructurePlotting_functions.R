
#############################

ggplotStructure=function(str,q){
	# ggplotStructure: 
	#	infile: input fastStructure table; assuming the second column has population labels
	#	qfile: table with .meanQ output; must have column headers
	require(ggplot2)
	pops=str[seq(1,nrow(str),2),2]
	inds=str[seq(1,nrow(str),2),1]
	struc5=q
	s5=stack(struc5)
	s5$indiv=rep(1:nrow(struc5))
#	s5$indiv=rep(inds)
	s5$pops=rep(pops)
	plo=ggplot(s5,aes(factor(indiv),values))+geom_bar(stat="identity",width=1,aes(fill=factor(ind)))+theme_bw()+scale_x_discrete(labels=pops)+xlab("pops")
	return(plo)
}

###############################
JSD.pair <- function(x, y){
	###Function to compute Shannon-Jensen Divergence
	###x and y are the frequencies for the same p categories
	u <- x/sum(x)
	v <- y/sum(y)
	m <- (u+v)/2
	if (all(u*v>0)){
		d <- (u*log(u/m)+v*log(v/m))/2
	} else {
		P1 <- u*log(u/m)
		P2 <- v*log(v/m)
		P1[is.nan(P1)] <- 0
		P2[is.nan(P2)] <- 0
		d <- (P1+P2)/2
	}
	return(sum(d))
}
##############################
matchPops=function(ga, gb, niter=3000) {
	### function to match population identifiers between fastStructure runs
	### based on permutations of column names and Shannon-Jensen divergences
	minsum=1000
	for (i in 1:niter) {
		names(gb)=sample(names(gb))
		sumjsd=0
		for (n in names(ga)) { 
			sumjsd=sumjsd+JSD.pair(ga[,n],gb[,n])
		}
		if (sumjsd<minsum) {
			minsum=sumjsd
			gbnames=names(gb)
		}
	}
	return(list("pops"=gbnames,"min.JSD"=minsum))
}
##############################
averageBest=function(likelihoods,top=25) {
	# matches populations assignments among best-likelihood runs,
	# averages assignemnt probabilities, returns averaged meanQ table
	bests=head(likes[order(likes[,2],decreasing=T),1],top)
	gs=read.table(paste(bests[1],".meanQ",sep=""))
	g1=gs
	print("top 1")
	for (b in 2:top) {
		gn=read.table(paste(bests[b],".meanQ",sep=""))
		names(gn)=matchPops(g1,gn)$pops
		print(paste("top",b))
		gs=gs+gn[,names(g1)]
	}
	return(gs/top)
}
###############################
