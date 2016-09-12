#' Dirichlet Process based High Dimensional Empirical Bayes Linear stimator
#'
#' This function is aimed to construct a linear classification rule based on multi-batch 
#' Dirichlet process empirical Bayes Classifier. Multi-batch 
#' Dirichlet process empirical Bayes is used to estimate 
#' mean difference to eliminate irrelavent features; then apply Naive Bayes classifier
#' 
#'
#' @param Train Matrix of training dataset.
#' @param Test  Matrix of testing dataset.
#' @param y Factor of class labels of training set
#' @param alpha Concentration parameter of Dirichlet process. 
#' @param sigma Standard deviation of Normal component in the base measure.
#' @param w Weight of zero component in the base measure.
#' @param T0 Upper bound of number of clusters. Default to be 10.
#' @param nfolds Number of batches we have to apply sparseVBDP function.
#' @param combine If combine is TRUE, we merge the component with small means (<0.5) into 
#' zero component. Default to be FALSE.
#' @param sparse A binary indicator. If sparse is TRUE, we use dis.CD.sparse function to
#' compute a sparse version of Empirical Bayes mean difference estimator; otherwise we use dis.CD function 
#' to compute Empirical Bayes mean estimator.
#' 
#' @return predicted class labels of observations in test dataset
#' @keywords VBDP
#' @export
#' @examples
#' #Leukemia data set classification
#' set.seed(100)
#' data(leukemia)
#' Train=as.matrix(leukemia[1:38,-1]);
#' Test=as.matrix(leukemia[-(1:38),-1]);
#' trainlabel=leukemia[1:38,1];
#' testlabel=leukemia[-(1:38),1];
#' table(DPclassifier(Train,Test,trainlabel,1,4,0.995),testlabel)


DPclassifier=function(Train,Test,y,alpha,sigma, w, T0=10,nfolds=10,combine=FALSE,sparse=TRUE){
  ###Train should be n times p matrix
  ###y should be a factor of 2 levels
  ###Test should be n0 times p matrix
  label=levels(y);
  if(length(label)!=2){stop("Please specify class labels for 2 groups!")}
  data1=Train[which(y==label[1]),];
  data2=Train[which(y==label[2]),]; 
  #compute sample size for each group
  n1=sum(y==label[1])
  n2=sum(y==label[2])
  mean1=apply(data1,2,mean)
  mean2=apply(data2,2,mean)
  var1=apply(data1,2,var)
  var2=apply(data2,2,var)
  S=sqrt(var1/n1+var2/n2)
  
  u1=t(t(data1)/S);
  u2=t(t(data2)/S);
  
  t=(mean1-mean2)/S
  
  
  results=multisparseVBDP(t,alpha,sigma, w, T0=T0,nfolds=nfolds,combine=combine)
  prior=results$prior;
  if(sparse==TRUE){vsparse=dis.CD.sparse(t,prior);}
  else{vsparse=dis.CD(t,prior);}
  a=vsparse/sqrt(sum(vsparse^2))
  a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;
  
  
  tdata=t(t(Test)/S);
  return(ifelse(tdata%*%a+a0>0,label[1],label[2]))
  
}










