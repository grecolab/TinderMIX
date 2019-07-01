#' This function construct confusion matrix between patient classes and the obtained clustering
#' @param classes is a vector of patient labels
#' @param clustering is a vector of clustering results
#' @param matrixRownames is a vector of names to assign as rownames of the confusion matrix
#' @param nCluster is the number of obtained clusters
#' @keywords multi-view clustering; confusion matrix
#' @return the confusion matrix 
#' @export
#' 
fisher_test <- function(classes, clustering, matrixRownames, nCluster){
  gene_clust <- clustering
  nClass <- length(table(classes))
  
  matrix_summary <-matrix(classes,ncol=1,nrow = length(classes));
  rownames(matrix_summary)<- matrixRownames;
  
  summary <- matrix_summary
  # for(i in 1:length(table(classes))){
  #   index <- which(matrix_summary[,1]==attr(table(classes)[i],which="name"));
  #   matrix_summary[index,1] <-i;
  # }
  
  ConfusionMatrix_gene <- matrix(0,nrow=nCluster,ncol=nClass)
  colnames(ConfusionMatrix_gene) = unique(classes)
  
  for(i in 1:nCluster){
    tab.i <- table(matrix_summary[which(gene_clust==i),1]);
    for(j in 1:length(tab.i)){
      #ConfusionMatrix_gene[i,as.integer(attr(tab.i[j],which="name"))] <- tab.i[j]
      ConfusionMatrix_gene[i,names(tab.i)[j]] <- tab.i[j]
      
    }
  }
  rownames(ConfusionMatrix_gene) <- paste("cluster",1:nCluster,sep="")
  #colnames(ConfusionMatrix_gene) <- paste("classi",1:nClass,sep="")
  
  CM_sup = ConfusionMatrix_gene
  colSums(CM_sup) -> elem_par_class
  
  fisherTest <- matrix(0,nrow=10,ncol=4)
  for(i in 1:nCluster){
    for(j in 1:nClass){
      cl1.test <- matrix(c(CM_sup[i,j],sum(CM_sup[i,-j]),elem_par_class[j],sum(elem_par_class[-j])),nrow=2,ncol=2,byrow=TRUE)
      fisherTest[i,j]<- fisher.test(cl1.test)$p.value
    }  
  }
  colnames(fisherTest) = colnames(CM_sup)
  rownames(fisherTest) = rownames(CM_sup)
  BFT = fisherTest<0.05
  
  return(list(CM_sup=CM_sup,fisherRes = fisherTest,BFT = BFT))
}