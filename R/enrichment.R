#'
#' This function perform enrichment for each gene label and compute a pathway wordcloud for every label
#'
#' @importFrom AnnotationDbi select
#' @importFrom gProfileR gprofiler
#' @importFrom wordcloud wordcloud
#' @importFrom gtools invalid
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#'
#' @param Mat matrix of gene labels
#' @param max.words max number of words in wordcluds
#' @param scale a vector of length 2 indicating the range of the size of the words.
#' @param random.order plot words in random order. If false, they will be plotted in decreasing frequency
#' @param min.freq words with frequency below min.freq will not be plotted
#' @param corrType string specifing the algorithm used for determining the significance threshold, one of gSCS, fdr, bonferroni. Default: fdr
#' @param type_enrich string specifying the enrichment type. Default = KEGG
#' @param org_enrich string specifying the organism. Default = rnorvegicus
#' @param pth numeric value specifyint the pvalue threshold. Default = 0.05
#' @param sig whether all or only statistically significant results should be returned
#' @param mis minimum size of functional category, smaller categories are excluded
#' @param only_annotated statistical domain size, one of "annotated", "known"
#' @return a list with the enriched pathways for each cluster of genes
#'
#' @export
#'

create_tic_tac_toe_wordcloud = function(Mat = res$Mat,max.words = 200,scale = c(0.8,2.5),random.order=FALSE,min.freq = 0,
                                        corrType = "fdr",type_enrich="REAC", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE){
  # library(wordcloud)
  # library(clusterProfiler)
  
  Enriched_list = list()
  for(i in c(1,4,7,2,5,8,3,6,9)){
    gi = Mat[,i]
    all_gi = names(gi[gi!=0])
    # gi_pos = names(gi[gi>0])
    # gi_neg = names(gi[gi<0])
    
    EP_all = compute_pathways(geneList = all_gi,corrType = corrType,type_enrich=type_enrich, org_enrich = org_enrich,pth = pth,sig = sig,mis = mis,only_annotated=only_annotated )
    # toRem = which(EP_all$Description %in% "Reactome")
    # if(length(toRem)>0) EP_all = EP_all[-toRem,]
    
    Enriched_list[[colnames(Mat)[i]]] = EP_all
    
  }
  
  
 par(mfrow = c(3,3), oma = c(0,0,0,0) + 0.1, mar = c(0,0,0,0) + 0.1)
  
  for(i in 1:length(Enriched_list)){
    EP_all = Enriched_list[[i]]
    if(nrow(EP_all)==0){
      plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    }else{
      #d = data.frame(word = substring(text = EP_all$Description,first = 1,last = 20), freq = log(EP_all$pValueAdj) * -2)
      d = data.frame(word = EP_all$Description, freq = log(EP_all$pValueAdj) * -2)
      
      # set.seed(1234)
      wordcloud(words = d$word, freq = d$freq, min.freq = min.freq,
                max.words=max.words, random.order=random.order, 
                colors=brewer.pal(8, "Dark2"), scale = scale)
      
    }
  }
  
  return(Enriched_list)
}

#'
#' This function create prototypes for a list of pathways
#'
#' @importFrom AnnotationDbi select
#' @importFrom gProfileR gprofiler
#' @importFrom gtools invalid
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#'
#' @param enrichedPath dataframe of enriched pathways coming from the compute_pathways function
#' @param annIDs vector with the pathways IDs for which the prototype will be computed
#' @param contour_res a list with the contours object returned in output by the create_contour function
#' @param nPerm number of permutation to run to compute the pvalue associated to pathway correaltion
#' @return a list with the prototypes of the pathways, their genes correlation and a pvalue
#'
#' @export
#'

create_pathway_prototypes = function(enrichedPath = enrichedPath, annIDs, contour_res = contour_res, nPerm = 100){
 # diss.cor <-1- abs(stats::cor(contour_res$GenesMap,method="pearson"))
  rownames(enrichedPath) = enrichedPath[,"annID"]
  RPGenes = contour_res$RPGenes
  names(RPGenes) = toupper(names(RPGenes))
  GenesMap = contour_res$GenesMap
  colnames(GenesMap) = toupper(colnames(GenesMap))
  
  Pathways_prot = list()
  
  for(i in 1:length(annIDs)){
    
    print(paste(i,"/", length(annIDs)))
    
    genes = unlist(strsplit(x = enrichedPath[annIDs[i],2],split = ","))
    pt = create_prot(RPGenes = RPGenes, genes=genes)
    cor_pt = mean(abs(stats::cor(GenesMap[,genes],method="pearson")))
    randomCor = c()
    pb = txtProgressBar(min = 1, max = nPerm, style = 3)
    for(permIdx in 1:nPerm){
      randomCor = c(randomCor, mean(abs(stats::cor(GenesMap[,sample(x = 1:ncol(GenesMap),size = length(genes))],method="pearson"))))
      setTxtProgressBar(pb,permIdx)
    }
    close(pb)
    pval_pt = 1 - sum(randomCor < cor_pt) / length(randomCor)
    protInfo = list(prototype = pt, correlation = cor_pt, pvalue = pval_pt)
    
    Pathways_prot[[enrichedPath[annIDs[i],"Description"]]] = protInfo
  }
  
  return(Pathways_prot)
}


create_prot = function(RPGenes,genes, nDosesInt = 50, nTimesInt = 50){
  mx = rep(0, nDosesInt)
  my = rep(0, nTimesInt)
  mz = matrix(0,nDosesInt,nTimesInt)
  for(j in 1:length(genes)){
    mz = mz + RPGenes[[genes[j]]][[3]]
  }
  mx = RPGenes[[genes[1]]][[1]]
  my = RPGenes[[genes[1]]][[2]]
  mz = mz/length(genes)
  L = list(mx, my, mz)
  return(L)
}


#'
#' This function perform enchment of the genes in each cluster
#'
#' @importFrom AnnotationDbi select
#' @importFrom gProfileR gprofiler
#' @importFrom gtools invalid
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#'
#' @param optimal_clustering vector of final clustering
#' @param corrType string specifing the algorithm used for determining the significance threshold, one of gSCS, fdr, bonferroni. Default: fdr
#' @param type_enrich string specifying the enrichment type. Default = KEGG
#' @param org_enrich string specifying the organism. Default = rnorvegicus
#' @param pth numeric value specifyint the pvalue threshold. Default = 0.05
#' @param sig whether all or only statistically significant results should be returned
#' @param mis minimum size of functional category, smaller categories are excluded
#' @param only_annotated statistical domain size, one of "annotated", "known"
#' @return a list with the enriched pathways for each cluster of genes
#'
#' @export
#'

compute_enrichment_for_clusters = function(optimal_clustering,corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE){

  nClust = length(unique(optimal_clustering))

  GList = list()
  for(i in 1:nClust){
    genelist = names(optimal_clustering)[optimal_clustering == i]
    GList[[i]] = cbind(genelist, 1:length(genelist))

  }

  GList = convert_genes(organism = org_enrich, GList=GList, annType = "SYMBOL")


  if(corrType == "none"){
    print("Nominal PValue")
    EnrichDatList = lapply(GList,enrich,type_enrich,org_enrich,pth,"bonferroni", sig = sig, mis = mis, only_annotated=only_annotated)
    for(i in 1:length(EnrichDatList)){
      ERi = EnrichDatList[[i]]
      ERi$pValueAdj = ERi$pValueAdj / length(ERi$pValueAdj)
      ERi$pValue = ERi$pValueAdj / length(ERi$pValue)
      EnrichDatList[[i]] = ERi
    }
  }else{
    EnrichDatList = lapply(GList,enrich,type_enrich,org_enrich,pth,corrType, sig = sig,  mis = mis, only_annotated=only_annotated)
  }

  return(EnrichDatList)
}

#'
#' This function perform enchment of a set of genes
#'
#' @importFrom AnnotationDbi select
#' @importFrom gProfileR gprofiler
#' @importFrom gtools invalid
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#'
#' @param geneList vector of gene identifiers
#' @param corrType string specifing the algorithm used for determining the significance threshold, one of gSCS, fdr, bonferroni. Default: fdr
#' @param type_enrich string specifying the enrichment type. Default = KEGG
#' @param org_enrich string specifying the organism. Default = rnorvegicus
#' @param pth numeric value specifyint the pvalue threshold. Default = 0.05
#' @param sig whether all or only statistically significant results should be returned
#' @param mis minimum size of functional category, smaller categories are excluded
#' @param only_annotated statistical domain size, one of "annotated", "known"
#' @return a list with the enriched pathways for each cluster of genes
#'
#' @export
#'
compute_pathways = function(geneList = rownames(res$Mat),corrType = "fdr",type_enrich="KEGG", org_enrich = "hsapiens",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE ){
  GList = list(matrix(geneList, ncol = 1))
  GList = convert_genes(organism = org_enrich, GList=GList, annType = "SYMBOL")
  
  if(corrType == "none"){
    print("Nominal PValue")
    EnrichDatList = lapply(GList,enrich,type_enrich,org_enrich,pth,"bonferroni", sig = sig, mis = mis, only_annotated=only_annotated)
    for(i in 1:length(EnrichDatList)){
      ERi = EnrichDatList[[i]]
      ERi$pValueAdj = ERi$pValueAdj / length(ERi$pValueAdj)
      ERi$pValue = ERi$pValueAdj / length(ERi$pValue)
      EnrichDatList[[i]] = ERi
    }
  }else{
    EnrichDatList = lapply(GList,enrich,type_enrich,org_enrich,pth,corrType, sig = sig,  mis = mis, only_annotated=only_annotated)
  }
  EnrichDatList = EnrichDatList[[1]]
  return(EnrichDatList)
  
}

#'
#' This function perform enchment of the genes in each cluster
#'
#' @importFrom AnnotationDbi select
#' @importFrom gProfileR gprofiler
#' @importFrom gtools invalid
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#'
#' @param x a dataframe with the gene names on the first column
#' @param adjust_method string specifing the algorithm used for determining the significance threshold, one of gSCS, fdr, bonferroni. Default: fdr
#' @param type string specifying the enrichment type. Default = KEGG
#' @param org string specifying the organism. Default = rnorvegicus
#' @param pval numeric value specifyint the pvalue threshold. Default = 0.05
#' @param sig whether all or only statistically significant results should be returned
#' @param mis minimum size of functional category, smaller categories are excluded
#' @param only_annotated statistical domain size, one of "annotated", "known"
#' @return a list with the enriched pathways for each cluster of genes
#'

enrich = function(x, type, org, pval, adjust_method,sig = FALSE, mis = 0, only_annotated = TRUE){
  if(only_annotated){
    domain_size = "annotated"
  }else{
    domain_size = "known"
  }

  print("before gprofiler")
  out = gProfileR::gprofiler(query = as.character(x[,1]),src_filter=type,organism=org,domain_size = domain_size,
                             max_p_value = pval, correction_method = adjust_method,significant = sig,min_isect_size = mis)

  #out = gProfileR::gprofiler(query = x[,1],src_filter=type,organism=org,max_p_value = pval,domain_size = "known")
  print(out)

  ##print("after gprofiler")
  out=out[,c("term.id","intersection","p.value","p.value","term.name")]
  colnames(out) = c( "annID","gID","pValue","pValueAdj","Description")

  print(out)

  if(nrow(out)>0)
    out$annID= gsub("KEGG:|REAC:","",out$annID)

  return(out)
}


#'
#' This function convert genes identifiers
#'
#' @importFrom AnnotationDbi select
#' @importFrom gtools invalid
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#'
#' @param organism a string specifying the organism under analysis
#' @param GList a list of genes identifier
#' @param annType string specifying the wanted gene identifier
#' @return a list with the converted genes identifiers
#'
convert_genes = function(organism = "hsapiens", GList, annType = "SYMBOL"){

  orgLibs <- list("hsapiens"=org.Hs.eg.db, "mouse"=org.Mm.eg.db, "rnorvegicus" = org.Rn.eg.db)
  orgDB <- orgLibs[[organism]]

  if(annType == "SYMBOL"){
    tmp = GList[[1]]
    tmp <- AnnotationDbi::select(orgDB, keys=as.character(tmp[,1]), columns="SYMBOL", keytype=annType)
    return(GList)
  }

  for(i in 1:length(GList)){
    M=GList[[i]]
    genes = M[,1]

    selectAnnDF <- AnnotationDbi::select(orgDB, keys=genes, columns="SYMBOL", keytype=annType)

    toRem = which(is.na(selectAnnDF[,2]))
    if(length(toRem)>0){
      selectAnnDF = selectAnnDF[-toRem,]
    }

    #remove eventual duplicates
    selectAnnDF = selectAnnDF[!duplicated(selectAnnDF$SYMBOL),]

    M = M[M[,1] %in% selectAnnDF[,1],]

    M = as.matrix(M)

    # macke sure they match by the key type
    matches = match(selectAnnDF[,1],M[,1])

    #all(M$ID[matches] == selectAnnDF[,1])
    M[matches,1] = selectAnnDF[,2]

    #update gene symbols
    GList[[i]] = M

  }

  return(GList)
}

#' This function create an excel file with the same format of the input need by the FunMappOne tool
#' @import xlsx
#'
#' @param optimal_clustering is a numeric vector with the clustering result for every gene
#' @param filePath is a string specifying the path of the xlsx file
#'
#' @export

write_xlsx_for_funmappone = function(optimal_clustering,filePath = "../contour_clustering/gene_clustering.xlsx"){
  genes1 = names(optimal_clustering)[optimal_clustering == 1]
  xlsx::write.xlsx(genes1, filePath, sheetName = "Cluster_1",row.names = FALSE)

  nClust = length(unique(optimal_clustering))

  for(i in 2:nClust){
    genesi = names(optimal_clustering)[optimal_clustering == i]
    xlsx::write.xlsx(genesi, filePath, sheetName = paste("Cluster_",i,sep=""),append = TRUE,row.names = FALSE)
  }

  grouping =  cbind(paste("Cluster_",1:nClust,sep=""), 1:nClust)
  colnames(grouping) = c("cluster", "group")
  xlsx::write.xlsx(grouping, filePath, sheetName = "grouping",append = TRUE, row.names = FALSE)

}
