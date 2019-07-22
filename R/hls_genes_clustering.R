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
#' It clusters the genes by their activity labels. It also estimate the ideal number of clusters
#' @importFrom stats cor lm hclust as.dist
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param GenesMap a matrix with the genes on the rows and the labels on the columns. The matrix contains 1, -1 or 0 if the gene has a label with an increasing or decreasing FC with respect to doses or no labels (0) 
#' @param nClust vector of putative numbers used to determine number of clusters. Default = c(5,10,25,50,75,100,125,150,175,200,250,300)
#' @param method string specifying the correlation method. Default = "pearson"
#' @param hls.method string specifying method for the hierarchical clustering. Default = "ward"
#' @return a list with the clustering results and statistics for the optimal number of clusters
#' \item{hls_res}{a list containing the clustering results, the clustering vector and the centers for each k value given in input}
#' \item{summaryMat}{a matrix with summary statistic for each k value given in input}
#' \item{clusterList}{a list with the clustering vectors}
#'
#' @export
#'
hls_labelbased_clustering = function(GenesMap, nClust = c(5,10,25,50,75,100,125,150,175,200,250,300), method="pearson", hls.method = "ward"){
  GenesMap = t(GenesMap)
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
#' @param optcl a vector with the optimal clustering 
#' @param mode string specifying the kind of prototype to be estimated. Possible values are mean, median or prot. mean and median compute the prototype as the mean of median values of the gene maps falling in the same cluster. Prot compute the prototype as the map of the most central gene in the cluster. default is mean
#' @return a list, named meanXYZ,  with the contour objects for each clustering prototype 
#'
#' @export
#'

create_prototypes = function(clust_res,contour_res, optcl, mode = "mean"){ #summaryMat
  
  if(mode %in% c("mean","median","prot")==FALSE){
    print("wrong mode parameter")
    return(NULL)
  }
  
  RPGenes = contour_res$RPGenes
  x = RPGenes[[1]][[1]]
  y = RPGenes[[1]][[2]]
  z = RPGenes[[1]][[3]]

  # pr = clust_res$hls_res[[which.max(summaryMat[,5])]]
  # optcl =pr$clusters
  meanXYZ = list()
  geneCluster = list()
  for(i in unique(optcl)){
    print(i)
    idxi = which(optcl == i)
    geneCluster[[i]] = names(optcl)[idxi]
    
    if(mode == "mean"){ # the cluster prototype is the mean value of the z-maps
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
    
    if(mode == "median"){ # the cluster prototype is the median value of the z-maps
      mx = RPGenes[[names(optcl)[1]]][[1]]
      my = RPGenes[[names(optcl)[1]]][[2]]
      mz = matrix(0, nrow = nrow(RPGenes[[names(optcl)[1]]][[3]]),ncol = ncol(RPGenes[[names(optcl)[1]]][[3]]))
    
      for(i in 1:nrow(mz)){
        for(j in 1:ncol(mz)){
          values_gij = c()
          for(gij in idxi){
            values_gij = c(values_gij,RPGenes[[names(optcl)[gij]]][[3]][i,j])
          }
          mz[i,j]=median(values_gij)
        }
      }
      meanXYZ[[i]] = list(mx, my, mz)
      
    }
    
    if(mode == "prot"){ # the cluster prototype is the most central gene in the cluster
      mx = RPGenes[[names(optcl)[1]]][[1]]
      my = RPGenes[[names(optcl)[1]]][[2]]
      
      mz = matrix(clust_res$hls_res[[which.max(clust_res$summaryMat[,5])]][[3]][i,], nrow(z),ncol(z))
      meanXYZ[[i]] = list(mx, my, mz)
      
    }

  }
  return(meanXYZ)
}

#' this function plots the clusters prototype
#' @importFrom plotly plot_ly
#'
#' @param meanXYZ a list the contour object for each prototype computed with the function create_prototypes
#' @param nR the number of rows to use in the plot
#' @param contour_size number specifying the distance between each contour line
#' 
#'
#' @export
#'

plot_clusters_prototypes = function(meanXYZ, nR = 2, nC = 5,
                                    activity_threshold = activity_threshold,
                                    BMD_resonse_threhold = BMD_resonse_threhold,
                                    nDoseInt = nDoseInt, nTimeInt = nTimeInt, 
                                    doseLabels = doseLabels, 
                                    timeLabels = timeLabels,
                                    oma = c(0,1,1,0),
                                    mar = rep(2,4)){
  
    label_leg = c()
    for(i in doseLabels){
      for(j in timeLabels){
        label_leg = c(label_leg, paste(i,j,sep="-"))
      }
    }

    #computing legend for each prototype
    labels = c()
    mainplots = list()     # this will be the list with the label of every prototype
    for(i in 1:length(meanXYZ)){
      print(i)
      immy = meanXYZ[[i]][[3]]
      coord = cbind(meanXYZ[[i]][[1]],meanXYZ[[i]][[2]])
      res2 = compute_BMD_IC50(immy,coord, i,
                              activity_threshold = activity_threshold,
                              BMD_resonse_threhold = BMD_resonse_threhold,
                              mode = mode,
                              nTimeInt = nTimeInt,nDoseInt=nDoseInt,
                              timeLabels = timeLabels,
                              doseLabels = doseLabels, toPlot = FALSE, tosave = FALSE)
      
      if(class(res2)== "numeric"){
        labels = rbind(labels, c(rep(0, nDoseInt*nTimeInt),0))
        if(res2 == 1){
          mainplots[[i]] = "No Response"
          
        }else{
          mainplots[[i]] = "No Dose Response"
        }
        
      }else{
        if(is.null(res2$label)){
          labels = rbind(labels, c(rep(0, nDoseInt*nTimeInt),0))
          mainplots[[i]] = "No Dose Response"
        }else{
          labels = rbind(labels, c(as.vector(res2$label$ttt_label),res2$verso))
          ttlab = paste(label_leg[as.vector(res2$label$ttt_label)!=0], collapse = " ")
          mainplots[[i]] = ttlab 
        }
        
      }
      
    }
    
    verso = labels[,ncol(labels)]
    labels = labels[,1:(ncol(labels)-1)]
    colnames(labels) = label_leg
    
    if(length(length(meanXYZ)) < (nR*nC)){
      m = matrix(1:(nR*nC), nrow = nR, ncol = nC, byrow = T)
      increased = FALSE
    }else{
      m = matrix(c(1:(nR*nC),rep((nR*nC)+1, nC)), nrow = nR + 1, ncol = nC, byrow = T)
      increased = TRUE
    }
    
    par(oma = oma + 0.1,
        mar = mar + 0.1)
    graphics::layout(mat = m)
    
    for(i in 1:length(meanXYZ)){
      print(i)
      immy = meanXYZ[[i]][[3]]
      coord = cbind(meanXYZ[[i]][[1]],meanXYZ[[i]][[2]])
      res2 = compute_BMD_IC50(immy,coord, strsplit(x = mainplots[[i]],split = " "),
                              activity_threshold = activity_threshold,
                              BMD_resonse_threhold = BMD_resonse_threhold,
                              mode = mode,
                              nTimeInt = nTimeInt,nDoseInt=nDoseInt,
                              timeLabels = timeLabels,
                              doseLabels = doseLabels, toPlot = TRUE, tosave = FALSE,addLegend = FALSE)
    }
    
    
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    graphics::legend(x="top",
                     c("Non responsive Area", "responsive Increasing","Responsive Decreasing",  "Dose-Response","IC50"),
                     col =c(NA, NA,NA,"gold", "red"),
                     lty = c(NA,NA,NA,1,1),
                     fill = c("darkblue","darkgreen","brown",NA,NA),
                     border = c("darkblue","darkgreen","brown",NA,NA),
                     lwd = c(NA,NA,NA,3,3),
                     xpd = NA, ncol = 1,box.lwd = 0,box.col = "white",bg = "white")
    
    
    # if(i<(nR*nC)){
    #   plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    #   graphics::legend(x="center",
    #                    c("Non responsive Area", "responsive Increasing","Responsive Decreasing",  "Dose-Response","IC50"),
    #                    col =c(NA, NA,NA,"gold", "red"),
    #                    lty = c(NA,NA,NA,1,1),
    #                    fill = c("darkblue","darkgreen","brown",NA,NA),
    #                    border = c("darkblue","darkgreen","brown",NA,NA),
    #                    lwd = c(NA,NA,NA,3,3),
    #                    xpd = NA, ncol = 1,box.lwd = 0,box.col = "white",bg = "white")
    #   
    # }else{
    #   while(i<(nR*nC)){
    #     plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    #     i = i+1
    #   }
    # }
    
    

    
}


# plot_clusters_prototypes_plotly = function(meanXYZ, nR = 2, contour_size=0.05){
# 
#   min_th = 100
#   max_th = -100
# 
#   for(i in 1:length(meanXYZ)){
#     mi = min(meanXYZ[[i]][[3]])
#     mix = max(meanXYZ[[i]][[3]])
#     if(mi<min_th){
#       min_th = mi
#     }
#     if(mix>max_th){
#       max_th = mix
#     }
#   }
# 
#   PL = list()
#   for(i in 1:length(meanXYZ)){
#     PL[[i]] <- plotly::plot_ly(x = meanXYZ[[i]][[1]],
#                                y = meanXYZ[[i]][[2]],
#                                z = t(meanXYZ[[i]][[3]]),
#                                type = "contour",
#                                autocontour = F,
#                                contours = list(
#                                  start = min_th,
#                                  end = max_th,
#                                  size = contour_size
#                                ))
#   }
#   p = plotly::subplot(PL,nrows = nR)
#   print(p)
# 
# }
