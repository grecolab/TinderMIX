library(TinderMIX)

data("FC_WY14643")
exp_data = fc_data#WY14643$exp_data
pheno_data = pdata#WY14643$pheno_data

print("Step 1: Compute Anova")

PvalMat = suppressMessages(compute_anova_dose_time(exp_data, pheno_data,dose_index = 2,time_point_index = 3))

print("Step 2: Identify responsive genes")
ItemsList = build_items_list(PvalMat)
# Responsive Genes are the gene that have a significant pvalue for dose, time and dose, time and dose * time
responsive_genes = unique(c(unlist(ItemsList$Dose),
                            unlist(ItemsList$Time),
                            unlist(ItemsList$`Dose:Time:DoseTime`),
                            unlist(ItemsList$`Dose:Time`)))

print("Step 3: Computing contour plot")
contour_res = suppressMessages(create_contour(exp_data, pheno_data, responsive_genes,dose_index = 2,time_point_index =3 ,gridSize = 50))
SST = contour_res$Stats
plot3d(toPlot = contour_res$RPGenes[["Pdk4"]],DF = contour_res$DFList[["Pdk4"]])
plot_clusters_prototypes(list(contour_res$RPGenes[["Nectin3"]]), nR=1)
plot_clusters_prototypes(list(contour_res$RPGenes[["Klf10"]]), nR=1)
plot_clusters_prototypes(list(contour_res$RPGenes[["Fam129a"]]), nR=1)
plot_clusters_prototypes(list(contour_res$RPGenes[["Dpp7"]]), nR=1)
plot_clusters_prototypes(list(contour_res$RPGenes[["Abcg8"]]), nR=1)

library(raster)
library(rasterVis)

rotate <- function(x) t(apply(x, 2, rev))
X = contour_res$RPGenes[["Klf10"]][[3]]
X = h(rotate(rotate(X)))
test.rast <- raster(X)
levelplot(test.rast)

write.table(contour_res$RPGenes[["Abcg8"]][[3]], file = "z_map_Abcg8.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
coord = cbind(contour_res$RPGenes[["Abcg8"]][[1]],contour_res$RPGenes[["Abcg8"]][[2]])
write.table(coord, file = "coord_Abcg8.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


write.table(contour_res$RPGenes[["Dpp7"]][[3]], file = "z_map_Dpp7.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
coord = cbind(contour_res$RPGenes[["Dpp7"]][[1]],contour_res$RPGenes[["Dpp7"]][[2]])
write.table(coord, file = "coord_Dpp7.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(contour_res$RPGenes[["Fam129a"]][[3]], file = "z_map_Fam129a.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
coord = cbind(contour_res$RPGenes[["Fam129a"]][[1]],contour_res$RPGenes[["Fam129a"]][[2]])
write.table(coord, file = "coord_Fam129a.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


write.table(contour_res$RPGenes[["Klf10"]][[3]], file = "z_map_Klf10.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
coord = cbind(contour_res$RPGenes[["Klf10"]][[1]],contour_res$RPGenes[["Klf10"]][[2]])
write.table(coord, file = "coord_Klf10.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

plot_clusters_prototypes(contour_res$RPGenes[1:10], nR=2)

# for(i in 1:length(contour_res$RPGenes)){
#   if()
# }

print("Step 4: Performing clustering")
hls_res = hls_genes_clustering(contour_res$GenesMap,  nClust = c(5,10,15,20,25), method="pearson", hls.method = "ward")

print("Step 5: Creating prototypes")
clpr = create_prototypes(clust_res = hls_res,summaryMat = hls_res$summaryMat,contour_res )

print("Step 6: print clustering prototype")
plot_clusters_prototypes(meanXYZ = clpr$meanXYZ, nR = 2)

GL = ItemsList[c("Dose","Time","Dose:Time:DoseTime","Dose:Time")]
classes = rep("D", length(responsive_genes))
names(classes) = responsive_genes
classes[names(classes) %in% GL$Time] = "T"
classes[names(classes) %in% GL$`Dose:Time`] = "DT"
classes[names(classes) %in% GL$`Dose:Time:DoseTime`] = "DTDT"

#clpr$optcl
res = fisher_test(classes = classes,clustering = hls_res$hls_res[[5]]$clusters,matrixRownames = names(classes) ,nCluster = 10)
View(res$BFT)

enrRes = compute_enrichment(clpr$optcl,corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE)
write_xlsx_for_funmappone(clpr$optcl,filePath = "../contour_clustering/gene_clustering.xlsx")
save.image("../contour_clustering/gene_clustering.RData")

CTD_C006253_pathways_20190616141646 <- read.csv("~/Downloads/CTD_C006253_pathways_20190616141646.csv")
enrRes[[1]]$Description %in% as.character(CTD_C006253_pathways_20190616141646$Pathway)
enrRes[[2]]$Description %in% as.character(CTD_C006253_pathways_20190616141646$Pathway)
enrRes[[3]]$Description %in% as.character(CTD_C006253_pathways_20190616141646$Pathway)
