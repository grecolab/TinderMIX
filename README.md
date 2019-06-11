# CUTENESS
The CUTENESS R package allows to analyse toxicogenomics gene expression dataset with multiple dose levels and time-points. It allows to identify the expression patterns with respect to both variables and to cluster genes accordingly. It also identify enriched pathways/go terms that are associated to each cluster.

# installing CUTENESS from GitHub

```R
library(devtools)
install_github("angy89/CUTENESS")
```

# Example of R Usage
```R
  library(CUTENESS)

  data("WY14643")
  exp_data = WY14643$exp_data
  pheno_data = WY14643$pheno_data

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
  plot_clusters_prototypes(list(contour_res$RPGenes[["Pdk4"]]), nR=1)

  print("Step 4: Performing clustering")
  hls_res = hls_genes_clustering(contour_res$GenesMap,  nClust = c(5,10,15,20,25), method="pearson", hls.method = "ward")

  print("Step 5: Creating prototypes")
  clpr = create_prototypes(clust_res = hls_res,summaryMat = hls_res$summaryMat,contour_res )

  print("Step 6: print clustering prototype")
  plot_clusters_prototypes(clpr$meanXYZ, nR = 2)

  enrRes = compute_enrichment(clpr$optcl,corrType = "fdr",type_enrich="KEGG", org_enrich = "rnorvegicus",pth = 0.05,sig = FALSE,mis = 0,only_annotated=FALSE)
  write_xlsx_for_funmappone(clpr$optcl,filePath = "../contour_clustering/gene_clustering.xlsx")
  ```

