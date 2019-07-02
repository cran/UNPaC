var.selection<-function(cluster_fun,cluster,k,x,p.adjust,alpha,clusterFirst=FALSE){
  selected_features=list()
  for (j in 2:k){
    if (clusterFirst==TRUE){cluster <- cluster_fun(x, j)$cluster}
  pvals<-apply(x,2,function(x) summary(aov(x~cluster))[[1]][1,5])
  selected_features[[j]]=which(p.adjust(pvals,p.adjust)<alpha)
  }
  return(selected_features)}
