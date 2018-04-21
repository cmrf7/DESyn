install.packages("limma")
install.packages("edgeR")
install.packages("fitdistrplus")
install.packages("matrixStats")
install.packages("locfit")
library(matrixStats)
library(edgeR)
library(fitdistrplus)

q.value=function(p, lambda=seq(0,0.95,0.05)) {
#
#This is Storey and Tibshirani's (PNAS, 2003) default method
#for determining the q-values.
#
#Code was originally obtained from John Storey's Web site.
#It has been edited slightly.
#
  m<-length(p) 
  pi0 <- rep(0,length(lambda))
  for(i in 1:length(lambda)) {
    pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])
  }
  spi0 <- smooth.spline(lambda,pi0,df=3)
  pi0 <- max(predict(spi0,x=1)$y,0)
  pi0 <- min(pi0,1)
  u <- order(p)
  v <- rank(p)
  qvalue <- pi0*m*p/v
  qvalue[u[m]] <- min(qvalue[u[m]],1)
  for(i in (m-1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
  }
  return(qvalue)
}


###"obs" is a matrix containing normalized counts from control and case groups; "gene.id" is a vector containing gene ids.

DESyn<-function(obs,gene.id){
repli<-dim(obs)[2]/2
cc<-dim(obs)[1]
###################################LRT test######################################################################
lrt.pvalue<-c(NA)
for (i in 1:cc){
if (all(obs[i,]==0)) {loglik0.overal<-0
}else {loglik0.overal<-fitdist(obs[i,],"nbinom",method="mle",control=list(maxit=5000))$loglik}
if (all(obs[i,1:repli]==0)) {loglik0.control<-0
}else {loglik0.control<-fitdist(obs[i,1:repli],"nbinom",method="mle",control=list(maxit=5000))$loglik}
if (all(obs[i,(repli+1):(repli*2)]==0)) {loglik0.case<-0
}else {loglik0.case<-fitdist(obs[i,(repli+1):(repli*2)],"nbinom",method="mle",control=list(maxit=5000))$loglik}
lrt.stat<--2*(loglik0.overal-loglik0.control-loglik0.case)
lrt.pvalue[i]<-1-pchisq(lrt.stat,df=2)
}
###################################Proposed########################################################################
proposed.pvalue<-c(NA)
proposed.stat<-matrix(c(NA),nrow=cc,ncol=(choose(repli*2,repli)/2))
for (j in 1:(choose(repli*2,repli)/2)){
grp<-rep(2,(repli*2))
grp[combn(repli*2,repli)[,j]]<-1
dm0<-DGEList(counts=obs,lib.size=rep(1,(repli*2)),norm.factors=rep(1,(repli*2)),group=grp)
disp.overal<-estimateDisp(dm0)$tagwise.dispersion
dm1<-DGEList(counts=obs[,grp==1],lib.size=rep(1,repli),norm.factors=rep(1,repli))
disp.control<-estimateDisp(dm1)$tagwise.dispersion
dm2<-DGEList(counts=obs[,grp==2],lib.size=rep(1,repli),norm.factors=rep(1,repli))
disp.case<-estimateDisp(dm2)$tagwise.dispersion
for (i in 1:cc){
if (all(obs[i,]==0)) {loglik.overal<-0
}else {loglik.overal<-fitdist(obs[i,],"nbinom",method="mle",fix.arg=list(size=1/disp.overal[i]),control=list(maxit=5000))$loglik}
if (all(obs[i,grp==1]==0)) {loglik.control<-0
}else {loglik.control<-fitdist(obs[i,grp==1],"nbinom",method="mle",fix.arg=list(size=1/disp.control[i]),control=list(maxit=5000))$loglik}
if (all(obs[i,grp==2]==0)) {loglik.case<-0
}else {loglik.case<-fitdist(obs[i,grp==2],"nbinom",method="mle",fix.arg=list(size=1/disp.case[i]),control=list(maxit=5000))$loglik}
proposed.stat[i,j]<--2*(loglik.overal-loglik.control-loglik.case)
}
}
null.stat<-proposed.stat[lrt.pvalue>=0.1,]
for (i in 1:cc){
proposed.pvalue[i]<-sum(null.stat>proposed.stat[i,1])/dim(null.stat)[1]/dim(null.stat)[2]
}

proposed.qvalue<-q.value(proposed.pvalue)
DE.gene<-gene.id[proposed.qvalue<0.05]
return(DE.gene)
}
