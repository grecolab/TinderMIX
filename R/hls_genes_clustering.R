#'
#' It clusters the contour plots by also estimating the ideal number of clusters
#' @importFrom stats cor lm hclust as.dist
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param GenesMap a matrix with the z-maps computed with the function create_contour
#' @param nClust vector of putative numbers used to determine number of clusters. Default = c(5,10,25,50,75,100,125,150,175,200,250,300)
#' @param method string specifying the correlation method. Default = "pearson"
#' @param hls.method string specifying method for the hierarchical clustering. Default = "ward"
#' @return a list with the clustering results and statistics for the optimal number of clusters
#' \item{hls_res}{a list containing the clustering results, the clustering vector and the centers for each k value given in input}
#' \item{summaryMat}{a matrix with summary statistic for each k value given in input}
#' \item{clusterList}{a list with the clustering vectors}
#'
#' @examples
#'   data("WY14643")
#'   exp_data = WY14643$exp_data
#'   pheno_data = WY14643$pheno_data
#'   PvalMat = compute_anova_dose_time(exp_data, pheno_data,dose_index = 2,time_point_index = 3)
#'   ItemsList = build_items_list(PvalMat)
#'   responsive_genes = unique(c(unlist(ItemsList$Dose),unlist(ItemsList$Time),unlist(ItemsList$`Dose:Time:DoseTime`),unlist(ItemsList$`Dose:Time`)))
#'   contour_res = create_contour(exp_data, pheno_data, responsive_genes,dose_index = 2,time_point_index =3 ,gridSize = 50)
#'   hls_res = hls_genes_clustering(contour_res$GenesMap,  nClust = c(5,10,15,20,25), method="pearson", hls.method = "ward")
#'
#' @export
#'

hls_genes_clustering = function(GenesMap, nClust = c(5,10,25,50,75,100,125,150,175,200,250,300), method="pearson", hls.method = "ward"){
  DB.diss <- 1-stats::cor(x=GenesMap,method = method)
  DD = stats::as.dist(DB.diss)

  hls_res  = list()
  index = 1
  summaryMat = c()
  clusterList = list()
  pb = txtProgressBar(min=0, max = length(nClust), style = 3)

  for(k in nClust){
    print(k)
    hls = stats::hclust(DD, method = hls.method)
    clusters = stats::cutree(hls,k = k)

    summaryMat = rbind(summaryMat,clustering_summary(DB = t(GenesMap),cluster = clusters))
    centers = findCenter(DB = t(GenesMap), clust_vector = clusters)

    hls_res[[index]] = list(hls = hls, clusters = clusters, centers = centers)
    clusterList[[index]] = clusters

    index = index + 1
    setTxtProgressBar(pb,index)
  }
  close(pb)

  return(list(hls_res = hls_res, summaryMat=summaryMat,clusterList=clusterList))
}

#'
#' This function evaluate the summary index for a clustering results by taking into account the within and between variability
#' @importFrom stats cor
#' @importFrom clv cls.scatt.diss.mx
#' @param DB your matrix dataset
#' @param cluster is a numeric vector of clustering results
#' @keywords clustering evaluation
#' @return a vector of evaluation indexes
#' @export

clustering_summary <- function(DB, cluster){
  diss.cor <-1- stats::cor(t(DB),method="pearson")
  toPrint <- matrix(0,nrow=1,ncol=5)
  NOBJS = dim(DB)[1]

  colnames(toPrint)<- c("#clusters","#singleton","Intra_avg","inter_complete","Index")
  tab <- table(cluster)
  center = findCenter(DB = DB,clust_vector = cluster)
  j = 1
  toPrint[j,1]<-length(tab) #number of clusters
  toPrint[j,2]<-length(which(tab==1)) #number of singleton

  clv::cls.scatt.diss.mx(diss.mx=diss.cor,clust=cluster) -> clust_eval
  clust_eval$intracls.average -> intra.avg
  clust_eval$intercls.complete -> inter.c

  toPrint[j,3] = round(1-mean(intra.avg),digits = 2)
  toPrint[j,4] = round(1-mean(inter.c),digits = 2)

  toPrint[j,5] =(1/3) * ( ((toPrint[j,3]+1)/2) + (1-(toPrint[j,4]+1)/2) +
                            (1-(toPrint[j,2]/length(table(cluster)))) ) # + toPrint[j,5])

  return(toPrint)
}

#'
#' This function select cluster prototypes. It select the prototype as the more correlated with the other elements in the cluster
#' @importFrom stats cor
#' @param DB your matrix dataset
#' @param clust_vector a numeric vector of clustering results
#' @keywords multi-view clustering; prototype; correlation
#' @return a matrix of prototypes
#' @export

findCenter <- function(DB,clust_vector){
  center <- NULL #matrice dei centroidi
  name_center <- c()
  elem_par_cluster <- table(clust_vector)
  for(i in names(elem_par_cluster)){
    if(elem_par_cluster[i]>1){
      clustElem <- DB[which(clust_vector==i),]
      matCor <- stats::cor(t(clustElem), method="pearson");
      bestCenter <- which.max(apply(matCor,1,FUN=function(riga){
        sum(riga)
      }))
      center <- rbind(center,clustElem[bestCenter,])
      name_center <- c(name_center,attr(bestCenter,which="names"))
    }else{
      clustElem <- DB[which(clust_vector==i),]
      center <- rbind(center,clustElem)
      name_center <- c(name_center,rownames(DB)[which(clust_vector==i)])
    }
  }
  rownames(center)<-name_center
  return(center)
}


#'
#' this function create the cluster prototypes as the mean values of the z maps of all the genes in the cluster
#' @param clust_res the clustering results object given in input by the function hls_genes_clustering
#' @param summaryMat the matrix with the summary statistics computed for the different k values
#' @param contour_res a list with the contours object returned in output by the create_contour function
#' @return a list with the contour objects for each clustering prototype and the vector of the optimal clustering
#' \item{meanXYZ}{a list the contour object for each prototype}
#' \item{optcl}{a vector with the optimal clustering}
#'
#' @examples
#'   data("WY14643")
#'   exp_data = WY14643$exp_data
#'   pheno_data = WY14643$pheno_data
#'   PvalMat = compute_anova_dose_time(exp_data, pheno_data,dose_index = 2,time_point_index = 3)
#'   ItemsList = build_items_list(PvalMat)
#'   responsive_genes = unique(c(unlist(ItemsList$Dose),unlist(ItemsList$Time),unlist(ItemsList$`Dose:Time:DoseTime`),unlist(ItemsList$`Dose:Time`)))
#'   contour_res = create_contour(exp_data, pheno_data, responsive_genes,dose_index = 2,time_point_index =3 ,gridSize = 50)
#'   hls_res = hls_genes_clustering(contour_res$GenesMap,  nClust = c(5,10,15,20,25), method="pearson", hls.method = "ward")
#'   clpr = create_prototypes(clust_res = hls_res,summaryMat = hls_res$summaryMat,contour_res )
#'
#' @export
#'

create_prototypes = function(clust_res,summaryMat,contour_res ){
  RPGenes = contour_res$RPGenes
  x = RPGenes[[1]][[1]]
  y = RPGenes[[1]][[2]]
  z = RPGenes[[1]][[3]]

  pr = clust_res$hls_res[[which.max(summaryMat[,5])]]

  optcl =pr$clusters
  meanXYZ = list()
  geneCluster = list()
  for(i in unique(optcl)){
    print(i)
    idxi = which(optcl == i)

    geneCluster[[i]] = names(optcl)[idxi]
    mx = rep(0, length(x))
    my = rep(0, length(y))
    mz = matrix(0,nrow(z),nrow(z))
    for(j in idxi){
      mx = mx + RPGenes[[names(optcl)[j]]][[1]]
      my = my + RPGenes[[names(optcl)[j]]][[2]]
      mz = mz + RPGenes[[names(optcl)[j]]][[3]]
    }
    mx = mx/length(idxi)
    my = my/length(idxi)
    mz = mz/length(idxi)
    meanXYZ[[i]] = list(mx, my, mz)
  }
  return(list(meanXYZ = meanXYZ, optcl = optcl))
}

#' this function plots the clusters prototype
#' @importFrom plotly plot_ly
#'
#' @param meanXYZ a list the contour object for each prototype computed with the function create_prototypes
#' @param nR the number of rows to use in the plot
#'
#' @examples
#'   data("WY14643")
#'   exp_data = WY14643$exp_data
#'   pheno_data = WY14643$pheno_data
#'   PvalMat = compute_anova_dose_time(exp_data, pheno_data,dose_index = 2,time_point_index = 3)
#'   ItemsList = build_items_list(PvalMat)
#'   responsive_genes = unique(c(unlist(ItemsList$Dose),unlist(ItemsList$Time),
#'                               unlist(ItemsList$`Dose:Time:DoseTime`),
#'                               unlist(ItemsList$`Dose:Time`)))
#'   contour_res = create_contour(exp_data, pheno_data, responsive_genes,\cr
#'                               dose_index = 2,time_point_index =3 ,gridSize = 50)
#'   hls_res = hls_genes_clustering(contour_res$GenesMap,  nClust = c(5,10,15,20,25),\cr
#'                                  method="pearson", hls.method = "ward")
#'
#'   clpr = create_prototypes(clust_res = hls_res,summaryMat = hls_res$summaryMat,contour_res )
#'   plot_clusters_prototypes(clpr$meanXYZ, nR = 2)
#'
#' @export
#'

plot_clusters_prototypes = function(meanXYZ, nR = 2){
  PL = list()
  for(i in 1:length(meanXYZ)){
    PL[[i]] <- plotly::plot_ly(x = meanXYZ[[i]][[1]], y = meanXYZ[[i]][[2]], z = t(meanXYZ[[i]][[3]]), type = "contour")
  }
  p = plotly::subplot(PL,nrows = nR)
  print(p)

}
