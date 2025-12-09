##' Spectrum Based Test
##'
##' This function employs the Spectrum based test statistic to infer
##' if covariance between two random vectors are zero or not.
##'
##' @param res1,res2  Two data matrice with the same sample size. The goal is to test
##' if the covariance between them are 0 or not.
##' @return Returns an list with \item{tstats}{ The test statistic.}
##' \item{pval}{The p-value computed based on the test statistic. The reference
##' distribution is weighted Chi-square. }
##' \item{weights}{The estimated weights for the weighted Chi-square of reference.}
##' @export

TestCCA=function(res1,res2){
  res1=scale(res1, center = TRUE, scale = FALSE)
  res2=scale(res2, center = TRUE, scale = FALSE)

  p1=dim(res1)[2];p2=dim(res2)[2];N=dim(res1)[1]

  mtn_cor = cor(cbind(res1,res2))[1:p1,(p1+1):(p1+p2)]

  W = cov(cbind(res1,res2))
  S11=W[1:p1,1:p1];S22=W[(p1+1):(p1+p2),(p1+1):(p1+p2)]
  S12=W[1:p1,(p1+1):(p1+p2)];S21=W[(p1+1):(p1+p2),1:p1]

  if (min(abs(eigen(S11)$values))<1.0e-8){
    ind_temp = which(abs(eigen(S11)$values)<1.0e-8)
    ev = eigen(S11)$values
    ev[ind_temp] = sqrt(log(max(p1,p2))/N)
    S11 = eigen(S11)$vectors%*%diag(ev)%*%t(eigen(S11)$vectors)
  }
  if (min(abs(eigen(S22)$values))<1.0e-8){
    ind_temp = which(abs(eigen(S22)$values)<1.0e-8)
    ev = eigen(S22)$values
    ev[ind_temp] = sqrt(log(max(p1,p2))/N)
    S22 = eigen(S22)$vectors%*%diag(ev)%*%t(eigen(S22)$vectors)
  }

  cca=svd(pracma::sqrtm(as.matrix(S11))$Binv%*%S12%*%solve(S22)%*%S21%*%pracma::sqrtm(as.matrix(S11))$Binv)$d[1:min(p1,p2)]
  if (cca[1]>0.99) {cca[which(cca>0.99)]=0.99}
  tstats=-N*log(prod(1-cca))

  #############
  chisqtail=function(lam,t) {
    c1=sum(lam);c2=sum(lam^2);c3=sum(lam^3);c4=sum(lam^4)
    s1=c3/c2^{1.5};s2=c4/c2^2
    a=1/(s1-sqrt(max(s1^2-s2,0)))
    del=round(s1*a^3-a^2,6)
    l=a^2-2*del

    ts=(t-c1)/sqrt(2*c2)
    ux=l+del
    sigx=sqrt(2)*a
    tail=pchisq(q=ts*sigx+ux,df=l,ncp=del,lower.tail=FALSE)

    return(tail)
  }

  res1=res1%*%pracma::sqrtm(as.matrix(S11))$Binv
  res2=res2%*%pracma::sqrtm(as.matrix(S22))$Binv

  dd=NULL
  for (i in 1:dim(res1)[2]) {

    for (j in 1:dim(res2)[2]){
      v=res1[,i]*res2[,j]
      dd=cbind(dd,v)
    }

  }

  Z=1/N*t(dd)%*%dd

  lam=eigen(Z)$values

  pval=chisqtail(lam,tstats)

  return(list(pval=pval,tstats=tstats,weights=lam))
}
