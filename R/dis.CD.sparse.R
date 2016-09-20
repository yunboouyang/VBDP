#' Compute sparse posterior mean based on discrete probability mass
#'
#' This function is aimed to compute posterior mean vector based on prior mass estimated by
#' sparseVBDP function. Additional sparsity is applied via thresholding posterior
#' probability being 0.
#' @param x Input n-dimensional noisy sparse Gaussian Sequence.
#' @param prior.mass Discrete n-dimensional prior mass vector. Each entry represents 1/n point mass at that entry.
#' @param h Error variance. Default to be 1. In real applications we assume it could be reliably estimated.
#' @param kappa Posterior probability threshold. If posterior probability of being 0 is larger
#' than kappa, then the resulted posterior mean estimate is 0. Default to be 0.5 corresponding to MAP estimator.
#' @return Posterior mean estimator of the high dimensional Gaussian sequence.
#' @keywords VBDP
#' @export
#' @examples
#' truemu=c(rep(0,900), rep(2,90),rep(10,10))
#' x=rnorm(1000, truemu)
#' results=sparseVBDP(x,1,6,0.01)
#' prior=results$prior
#' mu=dis.CD.sparse(x,prior)



dis.CD.sparse = function(x, prior.mass, h=1){
  # prior on mu is a discrete mixture over values provided in prior.mass
  # bandwith = h
  #check whether it is 0 or not
  n=length(prior.mass)
  A=as.data.frame(table(prior.mass));
  freq=A$Freq;
  uniq=as.numeric(levels(A$prior.mass))[A$prior.mass]
  
  if(sum(uniq==0)==0){
    tmp = outer(x, uniq, '-'); 
    tmp = exp(-tmp^2/(2*h));
    tmp=t(t(tmp)*freq)
    tmp = tmp/rowSums(tmp); 
    return(tmp %*% uniq)
  }else {
    location=which(uniq==0)
    tmp = outer(x, uniq, '-'); 
    tmp = exp(-tmp^2/(2*h));
    tmp=t(t(tmp)*freq)
    tmp = tmp/rowSums(tmp); 
    post=tmp %*% uniq;
    post[which(tmp[,location]>0.5)]=0;
    return(post)
  }
  
}

