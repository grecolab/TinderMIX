# TinderMIX

# Example of R Usage
```R
library(TinderMIX)
  
load("sample_data/cyclosporine AFC_data.rdata")
source("sample_data/set_parameters.R")

exp_data = dataFC$fc_data
pheno_data = dataFC$pdata
responsive_genes = rownames(exp_data)

contour_res = suppressMessages(create_contour(exp_data, pheno_data, 
                                              responsive_genes,
                                              dose_index = dose_index,
                                              time_point_index = time_point_index ,
                                              gridSize = gridSize,
                                              pvalFitting = pvalFitting,
                                              pvalFitting.adj.method = pvalFitting.adj.method,
                                              logScale = logScale, modelSelection = modelSelection))

res = run_all_BMD_IC50(contour_res = contour_res,
                       activity_threshold = activity_threshold,  
                       BMD_response_threshold = BMD_response_threshold, 
                       mode=mode, 
                       nDoseInt=nDoseInt, nTimeInt=nTimeInt, 
                       doseLabels = doseLabels, timeLabels = timeLabels,
                       tosave=FALSE, addLegend = FALSE, path = ".",
                       relGenes = contour_res$ggenes, toPlot = FALSE)

# example of map plot for Cd4 gene
geneName = "Cd4"
immy = contour_res$RPGenes[[geneName]][[3]]
coord = cbind(contour_res$RPGenes[[geneName]][[1]],contour_res$RPGenes[[geneName]][[2]])
res2 = compute_BMD_IC50(immy,coord, geneName,
                        activity_threshold = activity_threshold,
                        BMD_response_threshold = BMD_response_threshold,
                        mode = mode,
                        nTimeInt = nTimeInt,nDoseInt=nDoseInt,
                        timeLabels = timeLabels,
                        doseLabels = doseLabels,toPlot = T, addLegend = F)

# number of genes fitted for each model in the dose-responsive gene list
table(contour_res$Stats[rownames(res$MMA),"OptMod"])

# plot heatmap with number of genes for each label
plot_number_genes_labels(res,drugName = "cyclosporine A")

# plot radial plots with number of genes for the 12 combination of dose and time points
plot_cake_diagrams_time_dose_effect(res)


StatMat = cbind(res$MMA, contour_res$Stats[rownames(res$MMA),])
gi = rownames(StatMat)

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="rnorvegicus_gene_ensembl")

genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = gi, mart =ensembl)
toRem = which(genedesc$external_gene_name %in% names(table(genedesc$external_gene_name)[table(genedesc$external_gene_name)>1]))
toRem = toRem[seq_along(toRem) %%2 == 0]
if(length(toRem)>0)genedesc = genedesc[-toRem,]
rownames(genedesc) = genedesc$external_gene_name

labels = sapply(gi, FUN = function(gg)colnames(StatMat)[which(StatMat[gg,1:9]!=0)])
splitted_labels = do.call(rbind, strsplit(x = labels,split = "-"))
gene_sign = sapply(gi, FUN = function(gg)StatMat[gg,which(StatMat[gg,1:9]!=0)])

XX = cbind(gi, labels, splitted_labels, gene_sign)
XX  = cbind(XX, StatMat[rownames(XX), 10:20])

mmi = apply(XX[,6:8], 1, FUN = function(elem){
  if(elem[2]==1) dose_sign = "+" else dose_sign = "-"
  if(elem[1]==1) time_sign = "+" else time_sign = "-"
  if(elem[3]==1){
    dose_letter = "D"
    time_letter = "t"
  }
  if(elem[3]==2){
    dose_letter = "d"
    time_letter = "T"
  } 
  if(elem[3]==0){
    dose_letter = "D"
    time_letter = "T"
  }
  paste(dose_letter, dose_sign, time_letter, time_sign, sep="")
})

XX = cbind(mmi,genedesc[rownames(XX),2],XX)

Mat = res$MMA
load("updated_kegg_hierarhcy.rdata")
Enriched_list = list()
for(i in c(3,6,9,2,5,8,1,4,7)){
  gi = Mat[,i]
  all_gi = names(gi[gi!=0])
  EP_all = compute_pathways(geneList = all_gi,corrType = corrType,type_enrich=type_enrich, org_enrich = org_enrich,pth = pth,sig = sig,mis = mis,only_annotated=only_annotated )
  colnames(EP_all)[1] = "ID"
  EP_all = merge(EP_all, kegg_hierarchy,by.y = "ID") 
  Enriched_list[[colnames(Mat)[i]]] = EP_all
}  

PathList = list()
for(i in names(Enriched_list)){
  PathList[[i]] = as.character(Enriched_list[[i]][,"Pathway"])
}

intersect(PathList$`Intermediate-Late`, intersect(
  PathList$`Intermediate-Early`, intersect(PathList$`Sensitive-Middle`,
                                           intersect(PathList$`Sensitive-Late`, intersect(PathList$`Resilient-Late`,
                                                                                          intersect(PathList$`Resilient-Early`, intersect(PathList$`Resilient-Middle`, PathList$`Sensitive-Early`)))))
))

PathListTime = list("Early" = c(), "Middle" = c(),"Late" = c())
for(time in c("Early","Middle","Late")){
  idx_time = grep(pattern = time,ignore.case = TRUE,x = names(Enriched_list))
  for(i in idx_time){
    PathListTime[[time]] = c(PathListTime[[time]], as.character(Enriched_list[[i]][,"Pathway"]))
  }
  PathListTime[[time]] = unique(PathListTime[[time]]) 
}

library(gplots)
time_based_list = attr(which = "intersect",x = venn(PathListTime))
  ```

