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
#' @examples
#' data("FC_WY14643")
#' exp_data = fc_data
#' pheno_data = pdata
#'   PvalMat = compute_anova_dose_time(exp_data, pheno_data,dose_index = 2,time_point_index = 3)
#'   ItemsList = build_items_list(PvalMat)
#'   responsive_genes = unique(c(unlist(ItemsList$Dose),unlist(ItemsList$Time),unlist(ItemsList$`Dose:Time:DoseTime`),unlist(ItemsList$`Dose:Time`)))
#'   contour_res = create_contour(exp_data, pheno_data, responsive_genes,dose_index = 2,time_point_index =3 ,gridSize = 50)
#'   hls_res = hls_genes_clustering(contour_res$GenesMap,  nClust = c(5,10,15,20,25), method="pearson", hls.method = "ward")
#'   clpr = create_prototypes(clust_res = hls_res,summaryMat = hls_res$summaryMat,contour_res )
#'   enrRes = compute_enrichment(clpr$optcl,corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE)
#'
#' @export
#'

compute_enrichment = function(optimal_clustering,corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE){

  nClust = length(unique(optimal_clustering))

  GList = list()
  for(i in 1:nClust){
    genelist = names(optimal_clustering)[optimal_clustering == i]
    GList[[i]] = cbind(genelist, 1:length(genelist))

  }

  GList = convert_genes(organism = "Rat", GList=GList, annType = "SYMBOL")


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

  orgLibs <- list("Human"=org.Hs.eg.db, "Mouse"=org.Mm.eg.db, "Rat" = org.Rn.eg.db)
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
