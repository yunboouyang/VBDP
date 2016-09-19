# VBDP
R package to estimate sparse Gaussian sequence and construct linear classifier

##Maintainer
Yunbo Ouyang, Ph.D. candidate in Department of Statistics, University of Illinois at Urbana and Champaign

##Description
This package is aimed to estimate sparse Gaussian sequence and construct linear classifier based on this empirical Bayes estimator. 
To estimate Gaussian sequence, this package uses Dirichlet process mixture model to estimate the prior.


##Example
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
