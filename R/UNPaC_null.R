UNPaC_null_sig<- function(x, k,cluster.fun,nsim=100, rho=0.02, cov="glasso",center=center,scale=scale,WCSS=FALSE,kiter=20) {
## This version generates a new U every simulation
  ## Also does not center features prior to clustering
  out <- NA
  x.s=scale(x)
  if (cov=="glasso"){cur.scov=glasso::glasso(var(x.s), rho=rho)$w} else {cur.scov=var(x.s)}
  cur.iscov <- chol(cur.scov)
  cur.h1 <- apply(x, 2, function(x) find.h1(x,h.max=1000))
  x.var <- apply(x, 2, var)
  for (i in 1:nsim) {
    cur.ystar <- apply(x, 2, sample, replace=TRUE)
    xn <- scale(cur.ystar, scale=FALSE) +
      t(cur.h1*matrix(rnorm(nrow(x)*ncol(x)), ncol=nrow(x), nrow=ncol(x)))
    xn <- t((1/sqrt(1+cur.h1^2/x.var))*t(xn))
    xn <- t(colMeans(cur.ystar)+t(xn))
    U <- pnorm(matrix(rnorm(ncol(x)*row(x)), ncol = ncol(x)) %*% cur.iscov)
    S=sapply(1:ncol(x),function(i) quantile(xn[,i],U[,i]))
    xn.km <- cluster.fun(S, k)$cluster
   out[i] <- CI(S,xn.km)
          }
 return(out)
}

UNPaC_null_k<- function(x, k,cluster.fun,nsim=100, rho=0.02, cov="glasso",
                        center=center,scale=scale,WCSS=FALSE,kiter=20,d.power) {
  out <- matrix(NA,nrow=nsim,ncol=k)
  if (WCSS==TRUE){SSMat=matrix(NA,nrow=nsim,ncol=k)}
  x.s=scale(x)
  if (cov=="glasso"){cur.scov=glasso::glasso(var(x.s), rho=rho)$w} else {cur.scov=var(x.s)}
  cur.iscov <- chol(cur.scov)
  cur.h1 <- apply(x, 2, function(x) find.h1(x,h.max=1000))
  x.var <- apply(x, 2, var)
  for (i in 1:nsim) {
    cur.ystar <- apply(x, 2, sample, replace=TRUE)
    xn <- scale(cur.ystar, scale=FALSE) +
      t(cur.h1*matrix(rnorm(nrow(x)*ncol(x)), ncol=nrow(x), nrow=ncol(x)))
    xn <- t((1/sqrt(1+cur.h1^2/x.var))*t(xn))
    xn <- t(colMeans(cur.ystar)+t(xn))
    U <- pnorm(matrix(rnorm(ncol(x)*row(x)), ncol = ncol(x)) %*% cur.iscov)
    S=sapply(1:ncol(x),function(i) quantile(xn[,i],U[,i]))
    for (j in 1:k){
      xn.km <- cluster.fun(S, j)$cluster
      out[i,j] <- CI(S,xn.km)
      if (WCSS==TRUE){
        SSMat[i,j]=log(GapW.k(S,xn.km,d.power))
      }}
   }
  if (WCSS=="TRUE"){return(list(out=out,WCSS=SSMat))}else{return(list(out=out))}
}

 
CI<-function(x,cluster){
  wcss.feature <- numeric(ncol(x))
  for (k in unique(cluster)) {
    indices <- (cluster == k)
    if (sum(indices) > 1)
    wcss.feature <- wcss.feature + apply(scale(x[indices,], center = TRUE, scale = FALSE)^2, 2, sum)
  }
    wcss = sum(wcss.feature)
    Tss=sum(apply(scale(x, center = TRUE, scale = FALSE)^2, 2, sum))
    wcss/Tss
}
