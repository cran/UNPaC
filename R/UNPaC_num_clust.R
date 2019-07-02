#' Unimodal Non-Parametric Cluster (UNPaC) Test for Estimating Number of Clusters
#'
#' UNPaC for estimating the number of clusters Compares the cluster index (CI) from the original data to that
#' produced by clustering a simulated ortho-unimodal reference distribution generated using a Gaussian copula.
#' The CI is defined to be the sum of the within-cluster sum of squares about the cluster means divided by the total sum of squares.   
#' The number of clusters is chosen to maximize the difference between the data cluster index and the 
#' reference cluster indices, but additional rules are also implmented (See below). This method is similar to them method
#' described in Helgeson and Bair (2016) except a Gaussian copula approach is used to account for feature correlation and the rules for 
#' choosing the number of clusters are as described below.
#'
#' @param x a dataset with n observations (rows) and p features (columns) 
#' @param k maximum number of clusters considered. (default=10)
#' @param cluster.fun function used to cluster data. Function should return list containing a component "cluster." 
#' Examples include \code{\link[stats]{kmeans}} and  \code{\link[cluster]{pam}}. 
#' @param nsim a numeric value specifying the number of unimodal reference distributions used for testing (default=1000)
#' @param cov method used for approximating the covariance structure.  options include: "glasso" 
#' (See \code{\link{glasso}}), and  "est". (default = "glasso")
#' @param rho  a regularization parameter used in implementation of the graphical lasso. See documentation in  
#' @param scale should data be scaled such that each feature has variance equal to one prior to clustering 
#' (default=FALSE)
#' @param center should data be centered such that each feature has mean equal to zero prior to clustering 
#' (default=TRUE)
#' @param var_selection should dimension be reduced using feature filtering procedure? See description below. (default=FALSE)
#' @param p.adjust p-value adjustment method for additional feature filtering. See \code{\link[stats]{p.adjust}} 
#' for options. (default="fdr"). Not used if p.adjust="none."
#' @param gamma threshold for feature filtering procedure. See description below. Not used if var_selection=FALSE (default=0.10)
#' @param d.power Power in estimating the low of the within cluster dispersion for comparison to the Gap statistic. See  \code{\link[cluster]{clusGap}}. 
#'
#'
#'
#' @return
#' The function returns a list with the following components: 
#' \itemize{
#' \item{\code{BestK}}: {A matrix with 1 row and 4 columns named: "Max_CI","Max_CI_wi_1SE","Max_scaled_CI" and "Max_logWCSS_wi_1SE".
#' These correspond to the number of clusters, K, chosen by four different rules. "Max_CI choses K to maximize the difference in CI's between the true data and the
#' reference data. "Max_CI_wi_1SE" uses the "1-SE" criterion as in Tibshirani et al (2001), except for the CI.
#' "Max_scaled_CI" chooses K to maximize the difference in CIs from the observed and reference data scaled by the standard error of the reference data CIs.
#' "Max_logWCSS_wi_1SE" uses the Gap statistic and the "1-SE" criterion (Tibshirani et al, 2001) for choosing K.}
#' \item{\code{full_process}}: {A matrix containing the number of clusters, K, evaluated, the CI from the data, the average CI from the null 
#' distribution, the difference between the data CI and average null CI, the standard error for the difference in CIs, the log of the within cluster dispersion from the data, 
#' the average log of within cluster dispersion from the null data, The difference in within cluster dispersion (the Gap statistic), and the standard error for the Gap statistic.}
#' \item{\code{selected_features}}: {A vector of integers indicating the features retained by the feature filtering process.}}
#' 
#' 
#' 
#' @references
#' \itemize{
#'     \item Helgeson E and Bair E (2016). Non-Parametric Cluster Significance Testing with Reference to a Unimodal Null Distribution.
#'     arXiv preprint arXiv:1610.01424.
#'     \item Tibshirani, R., Walther, G. and Hastie, T. (2001). Estimating the number of data clusters via the Gap statistic. Journal of the Royal Statistical Society B, 63, 411-423.
#'     
#' }
#'
#' @details
#' There are two options for the covariance matrix used in generating the Gaussian 
#' copula: sample covariance estimation, \code{cov="est"}, which should be used if n>p, and the graphical lasso 
#' \code{cov="glasso"} which should be used if n<p. 
#' 
#' In high dimensional (n<p) settings a dimension reduction step can be implemented which selects features 
#' based on an F-test for difference in means across clusters. Features having a p-value less than a threshold 
#' \code{gamma} are retained. For additional feature filtering a p-value adjustment procedure (such as p.adjust="fdr") 
#' can be used. If no features are retained the resulting p-value for the cluster significance test is given as 1.
#'
#' @examples
#' 	 test1 <- matrix(rnorm(100*50), nrow=100, ncol=50)
#'   test1[1:30,1:50] <- rnorm(30*50, 2)
#'   test.edit<-scale(test1,center=TRUE,scale=FALSE)
#'   UNPaC_k<-UNPaC_num_clust(test.edit,k=5,kmeans,nsim=100,cov="est")
#' 	 
#' @export
#' @name UNPaC_num_clust
#' @author Erika S. Helgeson, David Vock, Eric Bair
#' 



UNPaC_num_clust<-function(x,k=10,cluster.fun, nsim=1000, cov="glasso", rho=0.02,
scale=FALSE,center=FALSE,var_selection=FALSE,p.adjust="none",gamma=.1,d.power=1){
  x<-scale(x,scale=scale,center=center)
  VS=NA
  if (var_selection==TRUE){
    VS=unique(unlist(var.selection(cluster.fun,k=k,x=x,p.adjust=p.adjust,gamma,clusterFirst=TRUE)))
    if (sum(!is.na(VS))==0){print("No selected features")
      return(list(BestK=rep(1,5),full_process=NA,SelectedFeatures=VS))
    }else{
    x=matrix(x[,VS],ncol=length(VS))
  }}
  
out<-matrix(nrow=k,ncol=9)
colnames(out)<-c("k","data_CI","ave_null_CI","CI_diff","CI_SE","data_logW","Null_logW","WCSS_Gap","logW_SE")

k_star=NA

nulldata<-UNPaC_null_k(x,nsim=nsim,cluster.fun=cluster.fun,cov=cov,k=k,rho=rho,
center=center,scale=scale,WCSS=TRUE,d.power=d.power)

for (i in 1:k){

out[i,1]<-i

cluster <- cluster.fun(x, i)$cluster

out[i,2] <- CI(x,cluster)

out[i,3]<-(1/nsim)*sum(nulldata[[1]][,i])

# difference
out[i,4]<-out[i,3]-out[i,2]

# sd
out[i,5]<-sd(nulldata[[1]][,i])*sqrt(1+1/nsim)

# GAP
out[i,6] <- log(GapW.k(x,cluster,d.power=d.power))

out[i,7]<-(1/nsim)*sum(nulldata[[2]][,i])

# difference
out[i,8]<-out[i,7]-out[i,6]

# sd
out[i,9]<-sd(nulldata[[2]][,i])*sqrt(1+1/nsim)
}

kmax=which.max( out[,4] )
CIMax1SE=Max1SE(out[,8],out[,9])
MaxScaled=which.max(c(0,out[2:k,4]/out[2:k,5]))
WCSSMax1SE=Max1SE(out[,8],out[,9])

BestK<-matrix(c(kmax,CIMax1SE,MaxScaled,WCSSMax1SE),ncol=4)
colnames(BestK)<-c("Max_CI","Max_CI_wi_1SE","Max_scaled_CI","Max_logWCSS_wi_1SE")

return(list(BestK=BestK,full_process=out,SelectedFeatures=VS))
}


GapW.k <- function(x,clus,d.power=d.power) {
  n <- nrow(x)
  ii <- seq_len(n)
  0.5 * sum(vapply(split(ii, clus), function(I) {
    xs <- x[I, , drop = FALSE]
    sum(dist(xs)^d.power/nrow(xs))
  }, 0))
}

Max1SE<-function(f, SE.f){
  K <- length(f)
  g.s <- f - SE.f
  if (any(mp <- f[-K] >= g.s[-1])) which.max(mp) else K
}
