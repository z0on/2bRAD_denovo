inf=commandArgs(T)

require(vegan)
ma = as.matrix(read.table(inf))
pp0=capscale(ma~1)
# eigenvectors
eigs=pp0$CA$eig

# selecting [non-outlier] samples within 0.05-0.95 quantile range of PC1
goodeigs=which(eigs<quantile(eigs,0.95) & eigs>quantile(eigs,0.05))

pp1=capscale(ma[goodeigs,goodeigs]~1)
pdf(file=paste(inf,"_pcoa_eigens.pdf",sep=""),height=4.8,width=8)
par(mfrow=c(1,2))
plot(pp1)
eigsn=pp1$CA$eig/sum(pp1$CA$eig)
plot(eigsn)
lines(bstick(nrow(ma)-1),col="red")
dev.off()
print(round(eigsn[1]/bstick(nrow(ma)-1)[1],2))

