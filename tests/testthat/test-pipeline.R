context("test-pipeline")

test_that("pipeline works", {
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
  
  print("Example of map analysis")
  geneName = "Fam129a"
  immy = contour_res$RPGenes[[geneName]][[3]]
  coord = cbind(contour_res$RPGenes[[geneName]][[1]],contour_res$RPGenes[[geneName]][[2]])
  res = compute_BMD_IC50(immy,coord, geneName,BMD_threshold = 0.58)
  
  print("Step 4: Performing clustering")
  hls_res = hls_genes_clustering(contour_res$GenesMap,  nClust = c(5,10,15,20,25), method="pearson", hls.method = "ward")

  print("Step 5: Creating prototypes")
  clpr = create_prototypes(clust_res = hls_res,summaryMat = hls_res$summaryMat,contour_res )

  print("Step 6: print clustering prototype")
  plot_clusters_prototypes(meanXYZ = clpr$meanXYZ, nR = 2)
  
  print("Step 7: identify cluster labels")
  create_tic_tac_toe(map = t(clpr$meanXYZ[[1]][[3]]))
  
  enrRes = compute_enrichment(clpr$optcl,corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE)
  #write_xlsx_for_funmappone(clpr$optcl,filePath = "../contour_clustering/gene_clustering2.xlsx")
  # save.image("../contour_clustering/gene_clustering.RData")

  # CTD_C006253_pathways_20190616141646 <- read.csv("~/Downloads/CTD_C006253_pathways_20190616141646.csv")
  # enrRes[[1]]$Description %in% as.character(CTD_C006253_pathways_20190616141646$Pathway)
  # enrRes[[2]]$Description %in% as.character(CTD_C006253_pathways_20190616141646$Pathway)
  # enrRes[[3]]$Description %in% as.character(CTD_C006253_pathways_20190616141646$Pathway)
  
})
