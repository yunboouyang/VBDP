#' Leukemia dataset to perform high dimensional classification
#'
#'  Leukemia data from high-density Affymetrix oligonucleotide arrays, 
#'  There are 7129 genes and 72 samples coming from two classes: 
#'  47 in class ALL (acute lymphocytic leukemia) and 25 in class 
#'  AML (acute mylogenous leukemia).
#'  
#'  
#'
#' @docType data
#'
#' @usage data(leukemia)
#'
#'
#' @keywords datasets
#'
#' 
#' @references GOLUB, T. R., et al. (1999). Molecular classifcation of cancer: Class discovery and class prediction by gene expression monitoring. Science 286: 531¨C537.
#' 
#' 
#' @source \href{http://portals.broadinstitute.org/cgi-bin/cancer/datasets.cgi}{Cancer Program Datasets}
#'
#' @examples
#' #Leukemia data set classification
#' set.seed(100)
#' data(leukemia)
#' Train=as.matrix(leukemia[1:38,-1]);
#' Test=as.matrix(leukemia[-(1:38),-1]);
#' trainlabel=leukemia[1:38,1];
#' testlabel=leukemia[-(1:38),1];
#' table(DPclassifier(Train,Test,trainlabel,1,4,0.995),testlabel)
#' 
"leukemia"