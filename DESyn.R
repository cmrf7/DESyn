DESyn<-function(data,group){
geneid<-data[,1]
obs<-data[,-1]
repli1<-sum(group==0)
repli2<-sum(group==1)
cc<-dim(obs)[1]
###################################LRT test######################################################################
lrt.pvalue<-c(NA)
for (i in 1:cc){
if (all(obs[i,]==0)) {loglik0.overal<-0
}else {loglik0.overal<-fitdist(obs[i,],"nbinom",method="mle",control=list(maxit=200000))$loglik}
if (all(obs[i,group==0]==0)) {loglik0.control<-0
}else {loglik0.control<-fitdist(obs[i,group==0],"nbinom",method="mle",control=list(maxit=200000))$loglik}
if (all(obs[i,group==1]==0)) {loglik0.case<-0
}else {loglik0.case<-fitdist(obs[i,group==1],"nbinom",method="mle",control=list(maxit=200000))$loglik}
lrt.stat<--2*(loglik0.overal-loglik0.control-loglik0.case)
lrt.pvalue[i]<-1-pchisq(lrt.stat,df=2)
}
###################################Proposed########################################################################
proposed.pvalue<-c(NA)
proposed.stat<-matrix(c(NA),nrow=cc,ncol=(choose(length(group),repli1)))
for (j in 1:(choose(length(group),repli1))){
grp<-rep(2,length(group))
grp[combn(length(group),repli1)[,j]]<-1
dm0<-DGEList(counts=obs,lib.size=rep(1,length(group)),norm.factors=rep(1,length(group)),group=grp)
disp.overal<-estimateDisp(dm0)$tagwise.dispersion
dm1<-DGEList(counts=obs[,grp==1],lib.size=rep(1,repli1),norm.factors=rep(1,repli1))
disp.control<-estimateDisp(dm1)$tagwise.dispersion
dm2<-DGEList(counts=obs[,grp==2],lib.size=rep(1,repli2),norm.factors=rep(1,repli2))
disp.case<-estimateDisp(dm2)$tagwise.dispersion
for (i in 1:cc){
if (all(obs[i,]==0)) {loglik.overal<-0
}else {loglik.overal<-fitdist(obs[i,],"nbinom",method="mle",fix.arg=list(size=1/disp.overal[i]),control=list(maxit=200000))$loglik}
if (all(obs[i,grp==1]==0)) {loglik.control<-0
}else {loglik.control<-fitdist(obs[i,grp==1],"nbinom",method="mle",fix.arg=list(size=1/disp.control[i]),control=list(maxit=200000))$loglik}
if (all(obs[i,grp==2]==0)) {loglik.case<-0
}else {loglik.case<-fitdist(obs[i,grp==2],"nbinom",method="mle",fix.arg=list(size=1/disp.case[i]),control=list(maxit=200000))$loglik}
proposed.stat[i,j]<--2*(loglik.overal-loglik.control-loglik.case)
}
}
null.stat<-proposed.stat[lrt.pvalue>=0.1,]
for (i in 1:cc){
pvalue[i]<-sum(null.stat>proposed.stat[i,1])/dim(null.stat)[1]/dim(null.stat)[2]
}

qvalue<-q.value(proposed.pvalue)

return(cbind(geneid,pvalue,qvalue))
}
