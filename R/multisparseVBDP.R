#' multi-batch Dirichlet Process based Empirical Bayes Estimator
#'
#' This function is aimed to estimate mean vector of sparse Gaussian sequence in the extremely
#' sparse case. Because of "rich gets richer" phenomenon, Simply using sparseVBDP function in extremely sparse case will
#' result in an estimator of all zero. To revise this estimator, we randomly
#' divide the whole vector into several folds and then apply sparseVBDP function. The 
#' revised estimator is used to construct high dimensional linear classification rule.
#'
#' @param x Input noisy sparse Gaussian Sequence.
#' @param alpha Concentration parameter of Dirichlet Process. 
#' @param sigma Standard deviation of Normal component in the base measure.
#' @param w Weight of zero component in the base measure.
#' @param T0 Upper bound of number of clusters. Default to be 10.
#' @param nfolds Number of batches we have to apply sparseVBDP function.
#' @param combine If combine is TRUE, we merge the component with small means (<0.5) into 
#' zero component. Default to be FALSE.
#' @return \emph{prior} Discrete prior vector containing each entry's estimated prior mean.
#' @return \emph{csize} The number of unique values in \emph{prior}.
#' @return \emph{prob} Probability of each element of unknown mean vector to be zero.
#' @keywords VBDP
#' @export
#' @examples
#' truemu=c(rep(0,900), rep(2,90),rep(10,10))
#' x=rnorm(1000, truemu)
#' results=multisparseVBDP(x,1,6,0.01)
#' prior=results$prior



multisparseVBDP=function(x,alpha,sigma, w, T0=10,nfolds=10,combine=FALSE){
  n=length(x);
  wholeprior=rep(0,n) #record whole prior 
  wholeprob=rep(0,n);
  foldid=sample(c(rep(seq(nfolds), times =n%/%nfolds),rep(1,n%%nfolds)))
  for(i in seq(nfolds)){
    which= foldid==i
    result=sparseVBDP(x[which],alpha,sigma, w, T0=T0,combine=combine)
    prior=result$prior;
    wholeprior[which]=prior
    prob=result$prob;
    wholeprob[which]=prob;
  }
  csize=length(unique(wholeprior));
  return(list(prior=wholeprior,csize=csize,prob=wholeprob));
}

