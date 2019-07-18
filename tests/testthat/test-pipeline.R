context("test-pipeline")

test_that("pipeline works", {
  library(TinderMIX)

  #data("FC_WY14643")
  load("data/FC_WY14643.RData")
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
  adj.pval = p.adjust(SST[,1],method = "bonferroni")
  SST = cbind(adj.pval,SST)
  ggenes = rownames(SST)[SST[,1]<0.05]
  
  print("Step 4: Run BMD_IC50 analysis on every gene")
  res = run_all_BMD_IC50(contour_res = contour_res,
                         activity_threshold = 0.58,  BMD_resonse_threhold = 0.95, 
                         nDoseInt=3, nTimeInt=3, 
                         doseLabels = c("Low","Mid","High"), timeLabels = c("High","Mid","Low"),
                         tosave=FALSE, path = ".",
                         relGenes = ggenes)
    
  
  print("Step 4: Performing clustering")
  hls_res = hls_genes_clustering(contour_res$GenesMap[,rownames(res$Mat)],  nClust = c(5,10,15,20,25), method="pearson", hls.method = "ward")

  print("Step 5: Creating prototypes")
  clpr = create_prototypes(clust_res = hls_res,summaryMat = hls_res$summaryMat,contour_res )

  print("Step 6: print clustering prototype")
  plot_clusters_prototypes(meanXYZ = clpr$meanXYZ, nR = 2)
  dim(clpr$meanXYZ)
  
  print("Step 7: identify cluster labels")
  label2DMap(map=t(clpr$meanXYZ[[2]][[3]]), th=0.95, nDoseInt=3, nTimeInt=3, doseLabels = c("Low","Mid","High"), timeLabels = c("High","Mid","Low"))
    
  enrRes = compute_enrichment(clpr$optcl,corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE)
  #write_xlsx_for_funmappone(clpr$optcl,filePath = "../contour_clustering/gene_clustering2.xlsx")
  # save.image("../contour_clustering/gene_clustering.RData")


})
