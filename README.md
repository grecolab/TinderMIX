# TinderMIX

TinderMIX is a new computational framework for dose- and time- dependent gene expression analysis which aims to identify groups of genes that show a dynamic dose-response behaviour.

# Example of R Usage
```R
library(TinderMIX)

pheno = read.delim("sample_data/WY-14643_pheno.txt")
exp_data = read.delim("sample_data/WY-14643_exp.txt")
drugName = "WY-14643"
fc_data = TinderMIX::compute_fc(exp_data = exp_data,pheno_data = pheno,dose_index = 2,time_index = 4)

responsive_genes = rownames(fc_data$fc_data)[1:100]
exp_data = fc_data$fc_data
pheno_data = fc_data$pdata

source("sample_data/set_parameters.R")

contour_res = suppressMessages(TinderMIX::create_contour(exp_data, pheno_data, 
		responsive_genes,
		dose_index = dose_index,
		time_point_index = time_point_index ,
		gridSize = gridSize,
		pvalFitting = pvalFitting,
		pvalFitting.adj.method = pvalFitting.adj.method,
		logScale = T, modelSelection = 1:3))

print("Step 4: Run BMD_IC50 analysis on every gene")
DDRGene = TinderMIX::run_all_BMD_IC50(contour_res = contour_res,
	 activity_threshold = activity_threshold,  
	 BMD_response_threshold = BMD_response_threshold, 
	 mode=mode, 
	 nDoseInt=nDoseInt, nTimeInt=nTimeInt, 
	 doseLabels = doseLabels, timeLabels = timeLabels,
	 tosave=FALSE, addLegend = FALSE, path = ".",
	 relGenes = contour_res$ggenes, toPlot = FALSE)

plot_dynamic_dose_responsive_map(contour_res, geneName="Gad1",activity_threshold,
																						BMD_response_threshold,mode,nTimeInt,nDoseInt,
																						timeLabels,doseLabels)

# plot heatmap with number of genes for each label
# source("R/plotting_functions.R")
gene_plot = TinderMIX::plot_number_genes_labels(res = DDRGene,drugName = drugName,timeLabels = timeLabels,doseLabels = doseLabels)
plot(gene_plot)

# plot radial plots with number of genes for the 12 combination of dose and time points
label_plot = TinderMIX::plot_cake_diagrams_time_dose_effect(res = DDRGene, timeLabels = timeLabels,doseLabels = doseLabels)

res = TinderMIX::create_gene_table(DDRGene,contour_res, nTimeInt,nDoseInt,biomart_dataset = "rnorvegicus_gene_ensembl")
  ```

