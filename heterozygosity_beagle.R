### if your input is ANGSD beagle genotype likelihoods (-doGlf 2) -> ...beagle.gz
fin = commandArgs(T)
glf <- read.table(fin, header=TRUE)[,-c(1:3)]
glf <- array(c(t(as.matrix(glf))), c(3, ncol(glf)/3, nrow(glf)))

# ### if your input is TEXT beagle genotype probs (-doGeno 8)
# glf <- read.table("myglf.geno.gz", header=FALSE)[,-c(1:2)]
# glf <- array(c(t(as.matrix(glf))), c(3, ncol(glf)/3, nrow(glf)))

# glf[,2,1:10] # CHECKME: geno likes across sfs bins (rows) and first 10 sites (cols) for sample 2

### in all cases, the SAF for a given site is just the genotype likelihoods ...
### from Nielsen et al 2012 PLoS One, we have SAF[0] = GL[0/0]/choose(2,0), SAF[1] = 2*GL[0/1]/choose(2,1), SAF[2] = GL[1/1]/choose(2,0)
### ... and of course everything cancels except the GL[./.]

### so: we can just use the typical EM algorithm from realSFS: 

EMstep <- function (sfs, GL) rowMeans(prop.table(sfs * GL, 2))

SFS <- matrix(1/3,3,dim(glf)[2])

maxiter <- 200
tol <- 1e-8

for(sample in 1:dim(glf)[2])
{
  for (iter in 1:maxiter)
  {
    upd <- EMstep (SFS[,sample], glf[,sample,])
    if (sqrt(sum((upd - SFS[,sample])^2)) < tol)
      break;
    SFS[,sample] <- upd
  }
  if (iter == maxiter) warning("increase maximum number of iterations")
}
print(c(fin,round(summary(SFS[2,1:300]),3)),quote=F)
save(SFS,file=paste(fin,"_zygosity.RData",sep="")) # SFS[2,] is estimated heterozygosity per individual

