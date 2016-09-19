# VBDP
R package to estimate sparse Gaussian sequence and construct linear classifier

##Maintainer
Yunbo Ouyang, Ph.D. candidate in Department of Statistics, University of Illinois at Urbana and Champaign

##Description
This package is aimed to estimate sparse Gaussian sequence and construct linear classifier based on this empirical Bayes estimator. 
To estimate Gaussian sequence, this package uses Dirichlet process mixture model to estimate the prior and then calculate posterior mean or posterior mean with posterior probability thresholding.


##Example
###Estimate sparse Gaussian sequence using posterior mean estimator
```{r}
truemu=c(rep(0,900), rep(2,90),rep(10,10))
x=rnorm(1000, truemu)
results=sparseVBDP(x,1,6,0.01)
prior=results$prior
mu=dis.CD(x,prior)
```
###Estimate sparse Gaussian sequence using posterior mean estimator with thresholding
```{r}
truemu=c(rep(0,900), rep(2,90),rep(10,10))
x=rnorm(1000, truemu)
results=sparseVBDP(x,1,6,0.01)
prior=results$prior
mu=dis.CD.sparse(x,prior)
```


###High dimensional classification of Leukemia dataset
```{r}
set.seed(100)
data(leukemia)
Train=as.matrix(leukemia[1:38,-1]);
Test=as.matrix(leukemia[-(1:38),-1]);
trainlabel=leukemia[1:38,1];
testlabel=leukemia[-(1:38),1];
table(DPclassifier(Train,Test,trainlabel,1,4,0.995),testlabel)
```

##Reference
Author's Manuscript  
