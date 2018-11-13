inf=commandArgs(T)

if (length(inf)==0) {
	message("siteSelector: selects shared sites in two-column tab-delimited sites tables for ANGSD\noutputs files with original name an extension .sel\nusage: Rscript siteSelector.R sitefile1 sitefile2 sitefile3 ... sitefileN\n")
	stop("specify input files\n")
}

inp=c();sites=c()
message("reading input files...")
for ( i in 1:length(inf)) {
	message(inf[i])
	inp[[i]]=read.table(inf[i],sep="\t")
	row.names(inp[[i]])=paste(inp[[i]][,1],inp[[i]][,2],sep=".")
	if (i==1) { 
		sites=row.names(inp[[i]]) 
	} else {
		sites=intersect(sites,row.names(inp[[i]]))
	}
}
message(paste(length(sites),"sites retained"))
message("writing output: ", paste(paste(inf,collapse="_"),".sel",sep=""))
inp[[1]]=inp[[1]][sites,]
write.table(inp[[1]],file=paste(paste(inf,collapse="_"),".sel",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

	