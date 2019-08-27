context("test-pipeline")

test_that("pipeline works", {
  library(TinderMIX)
  
  data("FC_WY14643")
  exp_data = fc_data#WY14643$exp_data
  pheno_data = pdata#WY14643$pheno_data
  
  # parameters for phenodata
  dose_index = 2
  time_point_index = 3
  
  #parameter for fitting
  gridSize = 50
  pvalFitting = 0.05
  pvalFitting.adj.method = "none"
  
  # parameter BMD identification
  activity_threshold = 0.1
  BMD_response_threshold = 0.5
  mode = "most_left"
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

  
      ######  Example of how to plot gene map
  geneName = "Asf1a"
  immy = contour_res$RPGenes[[geneName]][[3]]
  coord = cbind(contour_res$RPGenes[[geneName]][[1]],contour_res$RPGenes[[geneName]][[2]])
  toPlot = TRUE
  activity_threshold = 0.1
  
  res2 = compute_BMD_IC50(immy,coord, geneName,
                          activity_threshold = activity_threshold,
                          BMD_response_threshold = BMD_response_threshold,
                          mode = mode,
                          nTimeInt = nTimeInt,nDoseInt=nDoseInt,
                          timeLabels = timeLabels,
                          doseLabels = doseLabels,toPlot = T)
  print(res2$verso)
  
  print("Step 4: Run BMD_IC50 analysis on every gene")
  res = run_all_BMD_IC50(contour_res = contour_res,
                         activity_threshold = activity_threshold,  
                         BMD_response_threshold = BMD_response_threshold, 
                         mode=mode, 
                         nDoseInt=nDoseInt, nTimeInt=nTimeInt, 
                         doseLabels = doseLabels, timeLabels = timeLabels,
                         tosave=FALSE, addLegend = FALSE, path = ".",
                         relGenes = contour_res$ggenes, toPlot = FALSE)
  
  Mat = res$Mat
  Enriched_list = list()
  for(i in c(3,6,9,2,5,8,1,4,7)){
    gi = Mat[,i]
    all_gi = names(gi[gi!=0])
    # gi_pos = names(gi[gi>0])
    # gi_neg = names(gi[gi<0])

    EP_all = compute_pathways(geneList = all_gi,corrType = corrType,type_enrich=type_enrich, org_enrich = org_enrich,pth = pth,sig = sig,mis = mis,only_annotated=only_annotated )
    # toRem = which(EP_all$Description %in% "Reactome")
    # if(length(toRem)>0) EP_all = EP_all[-toRem,]

    Enriched_list[[colnames(Mat)[i]]] = EP_all

  }
  
  print("Step 5: compute pathways for every gene label and plot wordcloud")
  mat_enrich_res = create_tic_tac_toe_wordcloud(Enriched_list = Enriched_list, max.words = 400,scale = c(0.8,2.5),random.order=FALSE,min.freq = 0)
  
  library(reshape2)
  ELM = acast(mat_enrich_res$EL, word~label, value.var="pValueAdj")
  ELM[is.na(ELM)] = 1
  toRem = which(rownames(ELM) %in% "KEGG pathways")
  if(length(toRem)>0){
    ELM = ELM[-toRem,]
  }
  
  
  library(heatmaply)
  heatmaply(t(ELM), k_col = 2, k_row = 3, dendrogram = "none") %>% 
    layout(margin = list(l = 130, b = 40))
  
  mat_enrich_res$WCDF$freq = mat_enrich_res$WCDF$freq/10
  
  toRem = which(mat_enrich_res$WCDF$word %in% "KEGG pathways")
  if(length(toRem)>0){
    mat_enrich_res$WCDF = mat_enrich_res$WCDF[-toRem,]
  }
  library(ggwordcloud)
  set.seed(42)
  ggplot(
    mat_enrich_res$WCDF,
    aes(label = word, size = freq/10)
  ) +
    geom_text_wordcloud_area() +
    scale_size_area(max_size = 10) +
    theme_minimal() +
    facet_wrap(~label)
  
  print("Step 5 bis: compute pathways enriched from the dose-response genes")
      #####  compute pathways
  enrichedPath = compute_pathways(geneList = rownames(res$Mat),corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE )
      #####  create prototype for selected pathways (annIDs) and compute mean gene correlation and pvalue with permutation test
  PatProt = create_pathway_prototypes(enrichedPath = enrichedPath, annIDs = c("04010"),contour_res = contour_res,nPerm = 100)

      ##### example on how to plot pathway prototype
  immy = PatProt$`Carbohydrate digestion and absorption`$prototype[[3]]
  coord = cbind(PatProt$`Carbohydrate digestion and absorption`$prototype[[1]],PatProt$`Carbohydrate digestion and absorption`$prototype[[2]])
  res2 = compute_BMD_IC50(immy,coord, "Carbohydrate digestion and absorption",
                          activity_threshold = activity_threshold,
                          BMD_response_threshold = BMD_response_threshold,
                          mode = mode,
                          nTimeInt = nTimeInt,nDoseInt=nDoseInt,
                          timeLabels = timeLabels,
                          doseLabels = doseLabels)
  
  print("Step 6: Performing clustering based on gene maps")
  hls_res = hls_genes_clustering(GenesMap = contour_res$GenesMap[,rownames(res$Mat)],  nClust = nClust, method=method, hls.method = hls.method)

  print("Step 7: Identify optimal clustering based on contour plots, creating and plotting prototypes")
  pr = hls_res$hls_res[[which.max(hls_res$summaryMat[,5])]]
  optcl =pr$clusters
  clpr = create_prototypes(clust_res = hls_res,contour_res = contour_res,optcl=optcl, mode = "mean")   # mode can be, mean, median, or prot
  labels = plot_clusters_prototypes(meanXYZ = clpr, nR = 3, nC = 4,
                                    activity_threshold = activity_threshold,
                                    BMD_response_threshold = BMD_response_threshold,
                                    nDoseInt = nDoseInt, nTimeInt = nTimeInt, 
                                    doseLabels = doseLabels, 
                                    timeLabels = timeLabels, mode = mode)
  
  print("Step 7bis: Identify optimal clustering based on gene labels, creating and plotting prototypes")
  hls_res_labels = hls_labelbased_clustering(GenesMap=res$Mat, nClust = nClust, method="pearson", hls.method = "ward")
  pr_labels = hls_res_labels$hls_res[[which.max(hls_res$summaryMat[,5])]]
  optcl_labels =pr_labels$clusters
  clpr_labels = create_prototypes(clust_res = hls_res_labels,contour_res = contour_res,optcl=optcl_labels, mode = "mean")
  
  labels = plot_clusters_prototypes(meanXYZ = clpr_labels, nR = 4, nC = 3,
                                      activity_threshold = activity_threshold,
                                      BMD_response_threshold = BMD_response_threshold,
                                      nDoseInt = nDoseInt, nTimeInt = nTimeInt, 
                                      doseLabels = doseLabels, 
                                      timeLabels = timeLabels, mode = mode)
    
  print("Step 8: perform enrichment for each clusters in funmappone style")
  enrRes = compute_enrichment_for_clusters(optimal_clustering = optcl,corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE)
  #write_xlsx_for_funmappone(clpr$optcl,filePath = "../contour_clustering/gene_clustering2.xlsx")
  # save.image("../contour_clustering/gene_clustering.RData")


})
