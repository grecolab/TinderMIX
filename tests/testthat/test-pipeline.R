context("test-pipeline")

test_that("pipeline works", {
  
  data("FC_WY14643")
  exp_data = fc_data#WY14643$exp_data
  pheno_data = pdata#WY14643$pheno_data
  
  # parameters for phenodata
  dose_index = 2
  time_point_index = 3
  pvalFitting.adj.method = "bonferroni"
  
  #parameter for fitting
  gridSize = 50
  pvalFitting = 0.05
  
  # parameter BMD identification
  activity_threshold = 0.25
  BMD_resonse_threhold = 0.5
  mode = "cumulative"
  timeLabels = c("Late","Middle","Early")
  doseLabels = c("Sensitive","Intermediate","Resilient")
  nTimeInt = 3
  nDoseInt = 3
  
  #parameter clusterings
  nClust = c(5,10,15,20,25)
  method="pearson"
  hls.method = "ward"
  
  #parameter enrichment
  corrType = "fdr"
  type_enrich="KEGG"
  org_enrich = "rnorvegicus"
  pth = 0.05
  sig = FALSE
  mis = 0
  only_annotated=FALSE
  
  print("Step 1: Compute Anova")
  PvalMat = suppressMessages(compute_anova_dose_time(exp_data, 
                                                     pheno_data,
                                                     dose_index = dose_index,
                                                     time_point_index = time_point_index))

  print("Step 2: Identify responsive genes")
  ItemsList = build_items_list(PvalMat)
  
  # Responsive Genes are the gene that have a significant pvalue for dose, time and dose, time and dose * time
  responsive_genes = unique(c(unlist(ItemsList$Dose),
                              unlist(ItemsList$Time),
                              unlist(ItemsList$`Dose:Time:DoseTime`),
                              unlist(ItemsList$`Dose:Time`)))

  print("Step 3: Computing contour plot")
  contour_res = suppressMessages(create_contour(exp_data, pheno_data, 
                                                responsive_genes,
                                                dose_index = dose_index,
                                                time_point_index =time_point_index ,
                                                gridSize = gridSize,
                                                pvalFitting = pvalFitting,
                                                pvalFitting.adj.method=pvalFitting.adj.method))

  
  geneName = "Hspbp1" #"Acmsd
  immy = contour_res$RPGenes[[geneName]][[3]]
  #immy2 = clpr[[1]][[3]]
  coord = cbind(contour_res$RPGenes[[geneName]][[1]],contour_res$RPGenes[[geneName]][[2]])
  plot_contour_plot(immy, coord,geneName)
    
  
  res2 = compute_BMD_IC50(immy,coord, geneName,
                          activity_threshold = activity_threshold,
                          BMD_resonse_threhold = BMD_resonse_threhold,
                          mode = mode,
                          nTimeInt = nTimeInt,nDoseInt=nDoseInt,
                          timeLabels = timeLabels,
                          doseLabels = doseLabels)
  
  
  print("Step 4: Run BMD_IC50 analysis on every gene")
  res = run_all_BMD_IC50(contour_res = contour_res,
                         activity_threshold = activity_threshold,  
                         BMD_resonse_threhold = BMD_resonse_threhold, 
                         mode=mode, 
                         nDoseInt=nDoseInt, nTimeInt=nTimeInt, 
                         doseLabels = doseLabels, timeLabels = timeLabels,
                         tosave=FALSE, addLegend = FALSE, path = ".",
                         relGenes = contour_res$ggenes, toPlot = FALSE)
  

  enrichedPath = compute_pathways(geneList = rownames(res$Mat),corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE )
  PatProt = create_pathway_prototypes(enrichedPath = enrichedPath, annIDs = c("00270","04152","01100"),contour_res = contour_res, mode = "mean", allPath = FALSE)
  
  immy = PatProt$`Metabolic pathways`$prototype[[3]]
  coord = cbind(PatProt$`Metabolic pathways`$prototype[[1]],PatProt$`Metabolic pathways`$prototype[[2]])
  
  res2 = compute_BMD_IC50(immy,coord, "Metabolic pathways",
                          activity_threshold = activity_threshold,
                          BMD_resonse_threhold = BMD_resonse_threhold,
                          mode = mode,
                          nTimeInt = nTimeInt,nDoseInt=nDoseInt,
                          timeLabels = timeLabels,
                          doseLabels = doseLabels)
  
  print("Step 4: Performing clustering based on gene maps")
  hls_res = hls_genes_clustering(GenesMap = contour_res$GenesMap[,rownames(res$Mat)],  nClust = nClust, method=method, hls.method = hls.method)

  print("Step 5: Identify optimal clustering and creating prototypes")
  pr = hls_res$hls_res[[which.max(hls_res$summaryMat[,5])]]
  optcl =pr$clusters
  # mode can be, mean, median, or prot
  clpr = create_prototypes(clust_res = hls_res,contour_res = contour_res,optcl=optcl, mode = "mean")

  labels = plot_clusters_prototypes(meanXYZ = clpr, nR = 2, nC = 5,
                                    activity_threshold = activity_threshold,
                                    BMD_resonse_threhold = BMD_resonse_threhold,
                                    nDoseInt = nDoseInt, nTimeInt = nTimeInt, 
                                    doseLabels = doseLabels, 
                                    timeLabels = timeLabels)
  
  # print("Step 4bis: performing clustering based on gene labels")
  hls_res_labels = hls_labelbased_clustering(GenesMap=res$Mat, nClust = nClust, method="pearson", hls.method = "ward")
  # 
  print("Step 5: Identify optimal clustering and creating prototypes")
  pr_labels = hls_res$hls_res[[which.max(hls_res$summaryMat[,5])]]
  optcl_labels =pr$clusters
  clpr_labels = create_prototypes(clust_res = hls_res_labels,contour_res = contour_res,optcl=optcl_labels, mode = "mean")
  
  
  print("Step 6: print clustering prototype and identify cluster labels")
  #plot_clusters_prototypes_plotly(meanXYZ=clpr$meanXYZ)
  
  labels = plot_clusters_prototypes(meanXYZ = clpr_labels, nR = 4, nC = 4,
                                      activity_threshold = activity_threshold,
                                      BMD_resonse_threhold = BMD_resonse_threhold,
                                      nDoseInt = nDoseInt, nTimeInt = nTimeInt, 
                                      doseLabels = doseLabels, 
                                      timeLabels = timeLabels)
    

  enrRes = compute_enrichment_for_clusters(optimal_clustering = optcl,corrType = "fdr",type_enrich="KEGG", org_enrich = "hsapiens",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE)
  #write_xlsx_for_funmappone(clpr$optcl,filePath = "../contour_clustering/gene_clustering2.xlsx")
  # save.image("../contour_clustering/gene_clustering.RData")


})
