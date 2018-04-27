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
