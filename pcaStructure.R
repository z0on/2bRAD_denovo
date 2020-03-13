inf=commandArgs(T)

require(vegan)
ma = as.matrix(read.table(inf))
pp0=capscale(ma~1)
# eigenvectors
eigs=pp0$CA$eig

# selecting [non-outlier] samples within 0.05-0.95 quantile range of PC1
goodeigs=which(eigs<quantile(eigs,0.95) & eigs>quantile(eigs,0.05))

pp1=capscale(ma[goodeigs,goodeigs]~1)
# plot(pp1$CA$eig) 
# lines(bstick(pp1),col="red")
print(round(pp1$CA$eig[1]/bstick(pp1)[1],2))
