#' Sparse Dirichlet Process based Empirical Bayes Estimator using Variational Bayes Algorithm
#'
#' This function is aimed to estimate means of sparse Gaussian sequence. 
#' @param x Input noisy sparse Gaussian Sequence.
#' @param alpha Concentration parameter of Dirichlet Process. 
#' @param sigma Standard deviation of Normal component in the base measure.
#' @param w Weight of zero component in the base measure.
#' @param T0 Upper bound of number of clusters. Default to be 10.
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
#' results=sparseVBDP(x,1,6,0.01)
#' prior=results$prior



sparseVBDP=function(x,alpha, sigma, w, T0=10, combine=FALSE){
  n=length(x);
  lambda=1/sigma^2;
  phi0=matrix(rep(0,n*T0),nrow=n,ncol=T0);
  phi=matrix(rep(1/T0,n*T0),nrow=n,ncol=T0);
  while(max(abs(phi-phi0))>10^-3){
    phi0=phi;
    gamma1=1+colSums(phi)[1:(T0-1)];
    gamma2=alpha+rev(cumsum(rev(colSums(phi)))[1:(T0-1)]);
    tau1=(t(phi)%*%x)[1:T0,];
    tau2=lambda+colSums(phi);
    odds=log(w)-log(1-w)+log(sqrt(1/lambda*colSums(phi)+1))-tau1^2/(2*tau2);
    p=exp(odds)/(1+exp(odds))
    d1=c(digamma(gamma1)-digamma(gamma1+gamma2),0);
    d2=digamma(gamma2)-digamma(gamma1+gamma2);
    d3=d1+c(0,cumsum(d2));
    d=d3-0.5*(1-p)*((tau1/tau2)^2+1/tau2);
    S=outer(x,(1-p)*tau1/tau2,'*')+outer(rep(1,n),d,'*');
    E=exp(S);
    phi=E/rowSums(E);
  }
  mean=c(0,tau1/tau2);
  if(combine==TRUE){mean[which(abs(mean)<0.5)]=0;}
  
  newphi=cbind(phi%*%p,phi%*%diag(1-p));
  if(combine==TRUE){newphi[,1]=rowSums(newphi[,which(abs(mean)<0.5)]);}
  zeroprob=newphi[,1];   
  number=max.col(newphi);
  prior=mean[number];
  csize=length(unique(number));
  return(list(prior=prior,csize=csize,prob=zeroprob));
}





