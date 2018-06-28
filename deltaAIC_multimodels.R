# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf

weight.cutoff=1e-2 # lowest Akaike weight to report a model
npl=read.table("likes")

names(npl)=c("model","id","npara","ll")
#npl$model=as.factor(npl$model)

# finding minimum likelihood for each model
maxlike=c()
for (m in unique(npl$model)) {
	sub=subset(npl,model==m)
	maxlike=data.frame(rbind(maxlike,sub[sub$ll==max(sub$ll),]))
}
npara=maxlike$npara
likes=maxlike$ll

aic=2*npara-2*likes
delta.aic= aic-min(aic)
delta.aic[delta.aic>100]=100


# weights of evidence
exps=exp(-delta.aic/2)
wts=exps/sum(exps)
wts
maxlike$wts=wts
maxlike=maxlike[order(maxlike$wts,decreasing=T),]
maxlike=subset(maxlike,wts>weight.cutoff)
maxlike$model=factor(maxlike$model,levels=maxlike$model[order(maxlike$wts,decreasing=T)])
plot(wts~model,maxlike,las=2,log="y",xlab="",ylab="evidence weight")
maxlike
